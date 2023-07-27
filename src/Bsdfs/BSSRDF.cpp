#include "BSSRDF.hpp"
#include "Common/Texture.hpp"
#include "Common/Parallel.h"
#include "Fresnel.hpp"
#include "Sampler/Warp.hpp"

bool CatmullRomWeights(int size, const Float *nodes, Float x, int *offset,
                       Float *weights) {
    // Return _false_ if _x_ is out of bounds
    if (!(x >= nodes[0] && x <= nodes[size - 1])) return false;

    // Search for the interval _idx_ containing _x_
    int idx = FindInterval(size, [&](int i) { return nodes[i] <= x; });
    *offset = idx - 1;
    Float x0 = nodes[idx], x1 = nodes[idx + 1];

    // Compute the $t$ parameter and powers
    Float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;

    // Compute initial node weights $w_1$ and $w_2$
    weights[1] = 2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;

    // Compute first node weight $w_0$
    if (idx > 0) {
        Float w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        Float w0 = t3 - 2 * t2 + t;
        weights[0] = 0;
        weights[1] -= w0;
        weights[2] += w0;
    }

    // Compute last node weight $w_3$
    if (idx + 2 < size) {
        Float w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        Float w3 = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0;
    }
    return true;
}

Float SampleCatmullRom2D(int size1, int size2, const Float *nodes1,
                         const Float *nodes2, const Float *values,
                         const Float *cdf, Float alpha, Float u, Float *fval = nullptr,
                         Float *pdf = nullptr) {
    // Determine offset and coefficients for the _alpha_ parameter
    int offset;
    Float weights[4];
    if (!CatmullRomWeights(size1, nodes1, alpha, &offset, weights)) return 0;

    // Define a lambda function to interpolate table entries
    auto interpolate = [&](const Float *array, int idx) {
        Float value = 0;
        for (int i = 0; i < 4; ++i)
            if (weights[i] != 0)
                value += array[(offset + i) * size2 + idx] * weights[i];
        return value;
    };

    // Map _u_ to a spline interval by inverting the interpolated _cdf_
    Float maximum = interpolate(cdf, size2 - 1);
    u *= maximum;
    int idx =
            FindInterval(size2, [&](int i) { return interpolate(cdf, i) <= u; });

    // Look up node positions and interpolated function values
    Float f0 = interpolate(values, idx), f1 = interpolate(values, idx + 1);
    Float x0 = nodes2[idx], x1 = nodes2[idx + 1];
    Float width = x1 - x0;
    Float d0, d1;

    // Re-scale _u_ using the interpolated _cdf_
    u = (u - interpolate(cdf, idx)) / width;

    // Approximate derivatives using finite differences of the interpolant
    if (idx > 0)
        d0 = width * (f1 - interpolate(values, idx - 1)) /
             (x1 - nodes2[idx - 1]);
    else
        d0 = f1 - f0;
    if (idx + 2 < size2)
        d1 = width * (interpolate(values, idx + 2) - f0) /
             (nodes2[idx + 2] - x0);
    else
        d1 = f1 - f0;

    // Invert definite integral over spline segment and return solution

    // Set initial guess for $t$ by importance sampling a linear interpolant
    Float t;
    if (f0 != f1)
        t = (f0 - std::sqrt(std::max((Float)0, f0 * f0 + 2 * u * (f1 - f0)))) /
            (f0 - f1);
    else
        t = u / f0;
    Float a = 0, b = 1, Fhat, fhat;
    while (true) {
        // Fall back to a bisection step when _t_ is out of bounds
        if (!(t >= a && t <= b)) t = 0.5f * (a + b);

        // Evaluate target function and its derivative in Horner form
        Fhat = t * (f0 +
                    t * (.5f * d0 +
                         t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 +
                              t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
        fhat = f0 +
               t * (d0 +
                    t * (-2 * d0 - d1 + 3 * (f1 - f0) +
                         t * (d0 + d1 + 2 * (f0 - f1))));

        // Stop the iteration if converged
        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;

        // Update bisection bounds using updated _t_
        if (Fhat - u < 0)
            a = t;
        else
            b = t;

        // Perform a Newton step
        t -= (Fhat - u) / fhat;
    }

    // Return the sample position and function value
    if (fval) *fval = fhat;
    if (pdf) *pdf = fhat / maximum;
    return x0 + width * t;
}

Float IntegrateCatmullRom(int n, const Float *x, const Float *values,
                          Float *cdf) {
    Float sum = 0;
    cdf[0] = 0;
    for (int i = 0; i < n - 1; ++i) {
        // Look up $x_i$ and function values of spline segment _i_
        Float x0 = x[i], x1 = x[i + 1];
        Float f0 = values[i], f1 = values[i + 1];
        Float width = x1 - x0;

        // Approximate derivatives using finite differences
        Float d0, d1;
        if (i > 0)
            d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
        else
            d0 = f1 - f0;
        if (i + 2 < n)
            d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
        else
            d1 = f1 - f0;

        // Keep a running sum and build a cumulative distribution function
        sum += ((d0 - d1) * (1.f / 12.f) + (f0 + f1) * .5f) * width;
        cdf[i + 1] = sum;
    }
    return sum;
}



// BSSRDF Utility Functions
Float FresnelMoment1(Float eta) {
    Float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
            eta5 = eta4 * eta;
    if (eta < 1)
        return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 +
               2.49277f * eta4 - 0.68441f * eta5;
    else
        return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 -
               1.27198f * eta4 + 0.12746f * eta5;
}

Float FresnelMoment2(Float eta) {
    Float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
            eta5 = eta4 * eta;
    if (eta < 1) {
        return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 +
               0.07883f * eta4 + 0.04860f * eta5;
    } else {
        Float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
        return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 +
               458.843f * r_eta + 404.557f * eta - 189.519f * eta2 +
               54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
    }
}


Float SeparableBSSRDF::pdfSP(const SurfaceEvent & po, const Intersection & pi) const {
    const vec3 & ss = po.frame.bitTangent,ts  = po.frame.tangent,ns = po.frame.n;
    vec3 d = po.its->p - pi.p;
    vec3 dLocal(dot(ss, d), dot(ts, d), dot(ns, d));
    vec3 nLocal(dot(ss, pi.Ng), dot(ts, pi.Ng), dot(ns, pi.Ng));

    // Compute BSSRDF profile radius under projection along each axis
    Float rProj[3] = {std::sqrt(dLocal.y * dLocal.y + dLocal.z * dLocal.z),
                      std::sqrt(dLocal.z * dLocal.z + dLocal.x * dLocal.x),
                      std::sqrt(dLocal.x * dLocal.x + dLocal.y * dLocal.y)};

    // Return combined probability from all BSSRDF sampling strategies
    Float pdf = 0, axisProb[3] = {.25f, .25f, .5f};
    Float chProb = 1 / Float(3);
    for (int axis = 0; axis < 3; ++axis)
        for (int ch = 0; ch < 3; ++ch)
            pdf += pdfSr(ch, rProj[axis],po) * std::abs(nLocal[axis]) * chProb *
                   axisProb[axis];
    return pdf;
}

Spectrum SeparableBSSRDF::sampleS(const Scene & scene, Float u1, vec2 u2, const SurfaceEvent & po, Intersection * pi,
                                  Float * pdf) const {
    Spectrum s = sampleSp(scene, u1, u2, po, pi, pdf);
    pi->bsdf = new BSSRDFAdapater(this);
    return s;
}

Spectrum
SeparableBSSRDF::sampleSp(const Scene & scene, Float u1, const vec2 & u2, const SurfaceEvent & po, Intersection * pi,
                          Float * pdf) const {
    po.frame.tangent;
    vec3 ns = po.frame.n;
    vec3 ts = po.frame.tangent;
    vec3 ss = po.frame.bitTangent;
    vec3 vx, vy, vz;
    if ( u1 < .5f ) {
        vx = ss;
        vy = ts;
        vz = vec3(ns);
        u1 *= 2;
    } else if ( u1 < .75f ) {
        // Prepare for sampling rays with respect to _ss_
        vx = ts;
        vy = vec3(ns);
        vz = ss;
        u1 = ( u1 - .5f ) * 4;
    } else {
        // Prepare for sampling rays with respect to _ts_
        vx = vec3(ns);
        vy = ss;
        vz = ts;
        u1 = ( u1 - .75f ) * 4;
    }

    // Choose spectral channel for BSSRDF sampling
    int ch = clamp((int) ( u1 * 3 ), 0, 2);
    u1 = u1 * 3 - ch;

    Float r = sampleSr(ch, u2[0], po);
    if ( r < 0 )
        return Spectrum(0);
    Float phi = 2 * Constant::PI * u2[1];
    Float rMax = std::max(0.0015f,maxSr(ch,po));
    if ( r > rMax )
        return Spectrum(0);
    Float l = 2 * std::sqrt(sqr(rMax) - sqr(r));
    Intersection base;
    base.p = po.its->p + r * ( vx * std::cos(phi) + vy * std::sin(phi) ) -
             l * vz * 0.5f;
    vec3 pTarget = base.p + l * vz;

    /// may intersect more than one point
    struct IntersectionChain {
        std::optional < Intersection > si;
        IntersectionChain * next;
    };
    IntersectionChain * chain = new IntersectionChain();
    IntersectionChain * ptr = chain;
    std::optional < Intersection > its;
    Ray ray = base.spawnRay(pTarget);
    Float tHit = ray.farT;
    int hitCount = 0;
    std::stringstream fartS;
    Ray tempRay;
    while ( true ) {
        if(hitCount>0){

        }
        ptr->si = scene.intersect(ray);
        if ( ! ptr->si  )
            break;
        tHit -=ray.farT;
        if(tHit<0)
        {
           break;
        }
        base = ptr->si.value();
        if ( ptr->si->primitive == po.its->primitive ) {
            IntersectionChain * next = new IntersectionChain();
            ptr->next = next;
            ptr = next;
            hitCount ++;
            fartS<<ray.farT<<" ";
        }
        tempRay = ray;
//        ray.o = ptr->si->p;
      //  tHit -= ray.farT;
        ray = base.spawnRay(pTarget);

//        ray.farT = tHit;
    }

    // Randomly choose one of several intersections during BSSRDF sampling
    if ( hitCount == 0 )
        return Spectrum(0.0f);
    int selected = clamp((int) ( u1 * hitCount ), 0, hitCount - 1);
    int tempHitCount = hitCount;
    hitCount -= selected;
    auto s = fartS.str();
    while ( selected -- > 0 ) {
        IntersectionChain *  next = chain->next;
        delete chain;
        chain = next;
    }
    * pi = chain->si.value();
    while ( hitCount -- > 0 ) {
        IntersectionChain *  next = chain->next;
        delete chain;
        chain = next;
    }
    * pdf = this->pdfSP(po, * pi) / tempHitCount;
    // Compute sample PDF and sreturn the spatial BSSRDF term $\Sp$
    return Sr(0.01f,po);
    return this->Sr(distance(pi->p, po.its->p), po);
}


Spectrum DisneyBSSRDF::Sr(Float r, const SurfaceEvent & event) const {
    Spectrum R = color->eval(event.its);
    Spectrum d = scatterDistance->eval(event.its);
    r = ( r < 0.000001f ) ? 0.000001f : r;
    constexpr auto EIGHT_PI = 4.0f * Constant::TWO_PI;
    return R * ( exp( Spectrum( -r ) / d ) + exp( Spectrum( -r ) / ( 3.0f * d ) )  / ( EIGHT_PI * d * r ) );
}

float burley_max_cdf_calc( float max_r_d ){
    return 0.25f * ( 4.0f - std::exp( -max_r_d ) - 3.0f * exp( -max_r_d / 3.0f ) );
}
static float burley_max_r_d = 0.2f;
static float burley_max_cdf = burley_max_cdf_calc( burley_max_r_d );
static float burley_inv_max_cdf = 1.0f / burley_max_cdf;


Float DisneyBSSRDF::sampleSr(int ch, Float u, const SurfaceEvent & event) const {
    Spectrum R = color->eval(event.its);
    Spectrum d = scatterDistance->eval(event.its);
    constexpr auto quater_cutoff = 0.25f;

    // importance sampling burley profile
    const auto ret = ( u < quater_cutoff ) ? -d[ch] * log( 4.0f * u ) : -3.0f * d[ch] * log( ( u - quater_cutoff ) * 1.3333f );

    static const float OneMinusEpsilon = 0x1.fffffep-1;

    if (u < .25f) {
        // Sample the first exponential
        u = std::min<Float>(u * 4, OneMinusEpsilon);  // renormalize to [0,1)
        return d[ch] * std::log(1 / (1 - u));
    } else {
        // Second exponenital
        u = std::min<Float>((u - .25f) / .75f, OneMinusEpsilon);  // normalize to [0,1)
        return 3 * d[ch] * std::log(1 / (1 - u));
    }

    // ignore all samples outside the sampling range
    return ( ret > burley_max_r_d * d[ch] ) ? -1.0f : ret;
}

Float DisneyBSSRDF::pdfSr(int ch, Float r, const SurfaceEvent & event) const {
    Spectrum d = scatterDistance->eval(event.its);

    if (r < 1e-6f) r = 1e-6f;  // Avoid singularity at r == 0.

    // Weight the two individual PDFs as per the sampling frequency in
    // Sample_Sr().
    return (.25f * std::exp(-r / d[ch]) / (2 * Constant::PI * d[ch] * r) +
            .75f * std::exp(-r / (3 * d[ch])) / (6 * Constant::PI * d[ch] * r));
    // Sr(ch,r) = ( 0.25f * exp( -r / d[ch] ) / ( TWO_PI * d[ch] * r ) + 0.75f * exp( -r / ( 3.0f * d[ch] ) ) / ( SIX_PI * d[ch] * r )
    constexpr auto EIGHT_PI = 4.0f * Constant::TWO_PI;
    r = ( r < 0.000001f ) ? 0.000001f : r;
    return ( exp( -r / d[ch] ) + exp( -r / ( 3.0f * d[ch] ) ) ) / ( EIGHT_PI * d[ch] * r ) * burley_inv_max_cdf;
}

Float DisneyBSSRDF::maxSr(int ch, const SurfaceEvent & event) const {
    return  sampleSr(ch,0.999f,event);
    Spectrum d = scatterDistance->eval(event.its);
    return burley_max_r_d * d[ch];

}

Spectrum TabelBSSRDF::Sr(Float d, const SurfaceEvent & event) const {
    auto sigma_a = sigmaA->eval(event.its),sigma_s = sigmaS->eval(event.its);
    Spectrum  sigma_t = sigma_a + sigma_s,rho;
    for (int c = 0; c < 3; ++c)
        rho[c] = sigma_t[c] != 0 ? (sigma_s[c] / sigma_t[c]) : 0;

    Spectrum Sr(0.f);
    for (int ch = 0; ch < 3; ++ch) {
        // Convert $r$ into unitless optical radius $r_{\roman{optical}}$
        Float rOptical = d * sigma_t[ch];

        // Compute spline weights to interpolate BSSRDF on channel _ch_
        int rhoOffset, radiusOffset;
        Float rhoWeights[4], radiusWeights[4];
        if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(),
                               rho[ch], &rhoOffset, rhoWeights) ||
            !CatmullRomWeights(table.nRadiusSamples, table.radiusSamples.get(),
                               rOptical, &radiusOffset, radiusWeights))
            continue;

        // Set BSSRDF value _Sr[ch]_ using tensor spline interpolation
        Float sr = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                Float weight = rhoWeights[i] * radiusWeights[j];
                if (weight != 0)
                    sr += weight *
                          table.EvalProfile(rhoOffset + i, radiusOffset + j);
            }
        }

        // Cancel marginal PDF factor from tabulated BSSRDF profile
        if (rOptical != 0) sr /= 2 * Constant::PI * rOptical;
        Sr[ch] = sr;
    }
    // Transform BSSRDF value into world space units
    Sr *= sigma_t * sigma_t;
    return clamp(Sr,0,1);
}

Float TabelBSSRDF::sampleSr(int ch, Float u, const SurfaceEvent & event) const {
    auto sigma_a = sigmaA->eval(event.its),sigma_s = sigmaS->eval(event.its);
    Spectrum  sigma_t = sigma_a + sigma_s,rho;
    for (int c = 0; c < 3; ++c)
        rho[c] = sigma_t[c] != 0 ? (sigma_s[c] / sigma_t[c]) : 0;
    if (sigma_t[ch] == 0) return -1;
    auto t = SampleCatmullRom2D(table.nRhoSamples, table.nRadiusSamples,
                              table.rhoSamples.get(), table.radiusSamples.get(),
                              table.profile.get(), table.profileCDF.get(),
                              rho[ch], u) /
           sigma_t[ch];
    return  t;
}

Float TabelBSSRDF::pdfSr(int ch, Float r, const SurfaceEvent & event) const {
    auto sigma_a = sigmaA->eval(event.its),sigma_s = sigmaS->eval(event.its);
    Spectrum  sigma_t = sigma_a + sigma_s,rho;
    for (int c = 0; c < 3; ++c)
        rho[c] = sigma_t[c] != 0 ? (sigma_s[c] / sigma_t[c]) : 0;

    Float rOptical = r * sigma_t[ch];

    // Compute spline weights to interpolate BSSRDF density on channel _ch_
    int rhoOffset, radiusOffset;
    Float rhoWeights[4], radiusWeights[4];
    if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(), rho[ch],
                           &rhoOffset, rhoWeights) ||
        !CatmullRomWeights(table.nRadiusSamples, table.radiusSamples.get(),
                           rOptical, &radiusOffset, radiusWeights))
        return 0.f;

    // Return BSSRDF profile density for channel _ch_
    Float sr = 0, rhoEff = 0;
    for (int i = 0; i < 4; ++i) {
        if (rhoWeights[i] == 0) continue;
        rhoEff += table.rhoEff[rhoOffset + i] * rhoWeights[i];
        for (int j = 0; j < 4; ++j) {
            if (radiusWeights[j] == 0) continue;
            sr += table.EvalProfile(rhoOffset + i, radiusOffset + j) *
                  rhoWeights[i] * radiusWeights[j];
        }
    }

    // Cancel marginal PDF factor from tabulated BSSRDF profile
    if (rOptical != 0) sr /= 2 * Constant::PI * rOptical;
    return std::max((Float)0, sr * sigma_t[ch] * sigma_t[ch] / rhoEff);
}

Float TabelBSSRDF::maxSr(int ch, const SurfaceEvent & event) const {
    return sampleSr(ch,0.99,event);
}

Spectrum SeparableBSSRDF::Sw(const vec3 & w) const {
        Float c = 1 - 2 * FresnelMoment1(1 / eta);
        return Spectrum(1 - Fresnel::dielectricReflectance(1/eta, CosTheta(w))) / (c * Constant::PI);

}

TabelBSSRDF::TabelBSSRDF(Float eta, std::shared_ptr<Texture<Spectrum>> sigma_a,
                         std::shared_ptr<Texture<Spectrum>> sigma_s, Float g)
: SeparableBSSRDF(eta), table(100,64),sigmaA(sigma_a),sigmaS(sigma_s){
        computeBeamDiffusionBSSRDF(g,eta,&table);
}


Float BeamDiffusionMultScattering(Float sigmaS, Float sigmaA, Float g, Float eta, Float r) {
    const int nSamples = 100;

    ///forward-scattering light will mostly continue in the same direction
    ///backward-scattering light will mostly scatter back and absorb.
    sigmaS = sigmaS * (1-g);
    Float sigmaT = sigmaA + sigmaS;
    ///reduced albedo
    Float rhop = sigmaS / sigmaT;

    ///non-classical diffusion coefficient
    Float Dg = (2 * sigmaA + sigmaS)/( 3 * sqr(sigmaS + sigmaA));

    ///effective transport coefficient. model the multiple scattering inside the medium
    Float sigmaTr = sqrt(sigmaA/Dg);
    Float fm1 = FresnelMoment1(eta),fm2 = FresnelMoment2(eta);
    Float ze = -2 * Dg * (1+3*fm2)/(1-2*fm1);

    Float cPhi = .25f * (1-2*fm1);
    Float cE = 0.5f * (1-3 * fm2);

    Float Ed = 0;
    for(int i =0; i<nSamples ; i++){
        ///real point source and virtual source depth
        Float zr = - std::log(1-(i+0.5f)/nSamples) / sigmaT;
        Float zv = -zr + 2 *ze;


        Float dr = sqrt(sqr(r) + sqr(zr));
        Float dv = sqrt(sqr(r) + sqr(zv));


        Float phiD = Constant::INV_FOUR_PI / Dg * (std::exp(-sigmaTr * dr) / dr -
                                     std::exp(-sigmaTr * dv) / dv);
        Float EDn = Constant::INV_FOUR_PI * (zr * (1 + sigmaTr * dr) *
                              std::exp(-sigmaTr * dr) / (dr * dr * dr) -
                              zv * (1 + sigmaTr * dv) *
                              std::exp(-sigmaTr * dv) / (dv * dv * dv));
        Float E = phiD * cPhi + EDn * cE;
        Float kappa = 1 - std::exp(-2 * sigmaT * (dr + zr));
        Ed += kappa * rhop * rhop * E;
    }
    return Ed / nSamples;
}

// Media Inline Functions
inline Float PhaseHG(Float cosTheta, Float g) {
    Float denom = 1 + g * g + 2 * g * cosTheta;
    return Constant::INV_FOUR_PI * (1 - g * g) / (denom * std::sqrt(denom));
}

Float BeamDiffusionSingleScattering(Float sigmaS, Float sigmaA, Float g, Float eta, Float r) {
    Float sigma_t = sigmaA + sigmaS, rho = sigmaS / sigma_t;
    Float tCrit = r * std::sqrt(eta * eta - 1);
    Float Ess = 0;
    const int nSamples = 100;
    for (int i = 0; i < nSamples; ++i) {
        // Evaluate single scattering integrand and add to _Ess_
        Float ti = tCrit - std::log(1 - (i + .5f) / nSamples) / sigma_t;

        // Determine length $d$ of connecting segment and $\cos\theta_\roman{o}$
        Float d = std::sqrt(r * r + ti * ti);
        Float cosThetaO = ti / d;

        auto t = (1 - Fresnel::dielectricReflectance(-cosThetaO, eta));
        // Add contribution of single scattering at depth $t$
        Ess += rho * std::exp(-sigma_t * (d + tCrit)) / (d * d) *
               PhaseHG(cosThetaO, g) * (1 - Fresnel::dielectricReflectance(eta,-cosThetaO)) *
               std::abs(cosThetaO);
    }
    return Ess / nSamples;
}

void computeBeamDiffusionBSSRDF(Float g, Float eta, BSSRDFTable * t) {
    t->radiusSamples[0] = 0;
    t->radiusSamples[1] = 2.5e-3f;
    for (int i = 2; i < t->nRadiusSamples; ++i)
        t->radiusSamples[i] = t->radiusSamples[i - 1] * 1.2f;

    // Choose albedo values of the diffusion profile discretization
    for (int i = 0; i < t->nRhoSamples; ++i)
        t->rhoSamples[i] =
                (1 - std::exp(-8 * i / (Float)(t->nRhoSamples - 1))) /
                (1 - std::exp(-8));
    parallel_for([&](int i) {
        // Compute the diffusion profile for the _i_th albedo sample

        // Compute scattering profile for chosen albedo $\rho$
        for (int j = 0; j < t->nRadiusSamples; ++j) {
            Float rho = t->rhoSamples[i], r = t->radiusSamples[j];
            auto temp = 2 * Constant::PI * r * (BeamDiffusionMultScattering(rho, 1 - rho, g, eta, r) +
                                                BeamDiffusionSingleScattering(rho, 1 - rho, g, eta, r));
            if( isnan(temp)){

            }
            t->profile[i * t->nRadiusSamples + j] = temp;

        }

        // Compute effective albedo $\rho_{\roman{eff}}$ and CDF for importance
        // sampling
        t->rhoEff[i] =
                IntegrateCatmullRom(t->nRadiusSamples, t->radiusSamples.get(),
                                    &t->profile[i * t->nRadiusSamples],
                                    &t->profileCDF[i * t->nRadiusSamples]);
    }, t->nRhoSamples);
}

BSSRDFTable::BSSRDFTable(int nRhoSamples, int nRadiusSamples): nRhoSamples(nRhoSamples),
                                                               nRadiusSamples(nRadiusSamples),
                                                               rhoSamples(new Float[nRhoSamples]),
                                                               radiusSamples(new Float[nRadiusSamples]),
                                                               profile(new Float[nRadiusSamples * nRhoSamples]),
                                                               rhoEff(new Float[nRhoSamples]),
                                                               profileCDF(new Float[nRadiusSamples * nRhoSamples]) {}


Float BSSRDFAdapater::Pdf(const SurfaceEvent & event) const {
    if ( event.wi.z <= 0.0f || event.wo.z <= 0.0f )
        return 0.0f;
    return Warp::squareToCosineHemispherePdf(event.wo);}

Spectrum BSSRDFAdapater::sampleF(SurfaceEvent & event, const vec2 & u) const {
    event.wi = Warp::squareToCosineHemisphere(u);
    event.pdf = Warp::squareToCosineHemispherePdf(event.wi);
    event.sampleType = BXDFType(BSDF_REFLECTION | BSDF_DIFFUSE);

    return f(event);
}



Spectrum BSSRDFAdapater::f(const SurfaceEvent & event) const {
   return bssrdf->Sw(event.wi);
}

Float BSSRDFAdapater::eta(const SurfaceEvent & event) const {
    return bssrdf->eta;
}

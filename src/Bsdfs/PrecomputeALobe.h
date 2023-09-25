#pragma  once
#include "sampler/InterpolatedDistribution1D.hpp"
#include "Common/math.hpp"


#include <vector>
#include <memory>


template<int N>
class GaussLegendre
{
    std::array<float, N> _points;
    std::array<float, N> _weights;

    // Note: Actual integration is performed in floating point precision,
    // but the roots and weights are determined in double precision
    // before rounding. This fixes cancellation badness for higher-degree
    // polynomials. This is only done once during precomputation and
    // should not be much of an issue (hence why this code is also
    // not very optimized).
    double legendre(double x, int n)
    {
        if (n == 0)
            return 1.0;
        if (n == 1)
            return x;

        double P0 = 1.0;
        double P1 = x;
        for (int i = 2; i <= n; ++i) {
            double Pi = ((2.0*i - 1.0)*x*P1 - (i - 1.0)*P0)/i;
            P0 = P1;
            P1 = Pi;
        }
        return P1;
    }

    double legendreDeriv(double x, int n)
    {
        return n/(x*x - 1.0)*(x*legendre(x, n) - legendre(x, n - 1));
    }

    double kthRoot(int k)
    {
        // Initial guess due to Francesco Tricomi
        // See http://math.stackexchange.com/questions/12160/roots-of-legendre-polynomial
        double x = std::cos(Constant::PI*(4.0*k - 1.0)/(4.0*N + 2.0))*
                   (1.0 - 1.0/(8.0*N*N) + 1.0/(8.0*N*N*N));

        // Newton-Raphson
        for (int i = 0; i < 100; ++i) {
            double f = legendre(x, N);
            x -= f/legendreDeriv(x, N);
            if (std::abs(f) < 1e-6)
                break;
        }

        return x;
    }

public:
    GaussLegendre()
    {
        for (int i = 0; i < N; ++i) {
            _points[i] = float(kthRoot(i + 1));
            _weights[i] = float(2.0/((1.0 - sqr(_points[i]))*sqr(legendreDeriv(_points[i], N))));
        }
    }

    template<typename ValueType, typename F>
    inline ValueType integrate(F f)
    {
        ValueType result = f(_points[0])*_weights[0];
        for (int i = 1; i < N; ++i)
            result += f(_points[i])*_weights[i];
        return result;
    }

    int numSamples() const
    {
        return N;
    }

    const std::array<float, N> &points() const
    {
        return _points;
    }

    const std::array<float, N> &weights() const
    {
        return _weights;
    }
};




// Utility class used by the RoughDielectricBcsdf for precomputing
// the azimuthal scattering functions
    class PrecomputedAzimuthalLobe {
    public:
        static const  int AzimuthalResolution = 64;

    private:
        std::unique_ptr<vec3[]> _table;
        std::unique_ptr<InterpolatedDistribution1D> _sampler;
        float beta;
    public:
        PrecomputedAzimuthalLobe(std::unique_ptr<vec3[]> table);

        void sample(float cosThetaD, float xi, float &phi, float &pdf) const ;

        vec3 eval(float phi, float cosThetaD) const {


            float u = (AzimuthalResolution - 1) * phi * Constant::INV_TWO_PI;
            float v = (AzimuthalResolution - 1) * cosThetaD;
            int x0 = clamp(int(u), 0, AzimuthalResolution - 2);
            int y0 = clamp(int(v), 0, AzimuthalResolution - 2);
            int x1 = x0 + 1;
            int y1 = y0 + 1;
            u = clamp(u - x0, 0.0f, 1.0f);
            v = clamp(v - y0, 0.0f, 1.0f);

            return (_table[x0 + y0 * AzimuthalResolution] * (1.0f - u) + _table[x1 + y0 * AzimuthalResolution] * u) *
                   (1.0f - v) +
                   (_table[x0 + y1 * AzimuthalResolution] * (1.0f - u) + _table[x1 + y1 * AzimuthalResolution] * u) * v;
        }

        float pdf(float phi, float cosThetaD) const {
            float u = (AzimuthalResolution - 1) * phi * Constant::INV_TWO_PI;
            float v = (AzimuthalResolution - 1) * cosThetaD;
            return _sampler->pdf(v, int(u)) * float(AzimuthalResolution * Constant::INV_TWO_PI);
        }

        float weight(float cosThetaD) const {
            float v = (AzimuthalResolution - 1) * cosThetaD;
            return _sampler->sum(v) * (Constant::TWO_PI / AzimuthalResolution);
        }
    };




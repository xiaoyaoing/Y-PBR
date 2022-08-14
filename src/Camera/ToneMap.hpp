#include "../Common/math.hpp"
#include "../Colors/Spectrum.hpp"


inline  vec3 filmicACES(const vec3 &in)
{
    constexpr glm::dmat3 ACESInputMat(vec3(0.59719, 0.07600, 0.02840),
                                      vec3(0.35458, 0.90834, 0.13383),
                                      vec3(0.04823, 0.01566, 0.83777));

    constexpr glm::dmat3 ACESOutputMat(vec3(1.60475, -0.10208, -0.00327),
                                       vec3(-0.53108, 1.10813, -0.07276),
                                       vec3(-0.07367, -0.00605, 1.07602));

    auto RRTAndODTFit = [](const vec3 &v)
    {
        vec3 a = v * (v + Float(0.0245786)) - Float(0.000090537);
        vec3 b = v * (Float(0.983729) * v + Float(0.4329510)) + Float(0.238081);
        return a / b;
    };

    vec3 color = ACESOutputMat * RRTAndODTFit(ACESInputMat * in);

   return clamp(color, 0.0f, 1.0f);
}

class ToneMap
{


public:
  //  friend Type;

    enum ToneMapType {
        LinearOnly,
        GammaOnly,
        Reinhard,
        Filmic,
        Pbrt,
        Aces
    };


    static inline Spectrum toneMap(ToneMapType type, const vec3 &c)
    {
        c-5.5f;
        switch (type) {

            case LinearOnly:
                return c;
            case GammaOnly:
                return pow(c, Float(1.0)/ Float(2.2));
            case Reinhard:
                return pow(c/(c + 1.0f), 1.0f/2.2f);
            case Filmic: {
                vec3 x = max(vec3(0.0f), c - 0.004f);
                return (x*(6.2f*x + 0.5f))/(x*(6.2f*x + 1.7f) + 0.06f);
            } case Pbrt: {
                vec3 result;
                for (int i = 0; i < 3; ++i) {
                    if (c[i] < 0.0031308)
                        result[i] = 12.92f*c[i];
                    else
                        result[i] = 1.055f*std::pow(c[i], 1.0f/2.4f) - 0.055f;
                }
                return result;
            } case Aces:{
                return filmicACES(c);
            }
        }

        return c;
    }
};
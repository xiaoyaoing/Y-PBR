//#include "ToneMap.hpp"
//
//Spectrum ToneMap::toneMap(ToneMap::ToneMapType type, const vec3 &c)
//    {
//
//        switch (type) {
//
//            case LinearOnly:
//                return c;
//            case GammaOnly:
//                return pow(c, Float(1.0)/ Float(2.2));
//            case Reinhard:
//                return pow(c/(c + 1.0f), 1.0f/2.2f);
//            case Filmic: {
//                vec3 x = max(vec3(0.0f), c - 0.004f);
//                return (x*(6.2f*x + 0.5f))/(x*(6.2f*x + 1.7f) + 0.06f);
//            } case Pbrt: {
//                vec3 result;
//                for (int i = 0; i < 3; ++i) {
//                    if (c[i] < 0.0031308)
//                        result[i] = 12.92f*c[i];
//                    else
//                        result[i] = 1.055f*std::pow(c[i], 1.0f/2.4f) - 0.055f;
//                }
//                return result;
//            } case Aces:{
//                return filmicACES(c);
//            }
//        }
//
//        return c;
//    }
//
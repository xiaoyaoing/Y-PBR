#include "Common/Texture.hpp"

template<class T>
class ConstantTexture : public  Texture<T>{
public:
    ConstantTexture(const T  value):value(value){
    }

    T Evaluate(const Intersection * si = nullptr) const override { return value; }
    T Evaluate(const vec2 & uv) const override {return value; }
    vec2 sample(TextureMapJacobian jacobian, const vec2 & uv,Float * pdf) const override {
        return vec2();
    }
    Float pdf(TextureMapJacobian jacobian, const vec2 & uv) const override {
        return 1;
    }

    void setScale(Float scale) override {
        value *= scale;
    }

protected:
    T value;
};
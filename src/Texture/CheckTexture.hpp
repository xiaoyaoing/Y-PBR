#include "Common/Texture.hpp"

template<class T>
class CheckTexture : public  Texture<T>{
public:
    CheckTexture(const  T & onColor,const T & offColor,int resU,int resV )
        :_onColor(onColor),_offColor(offColor),_resU(resU),_resV(resV){
    }


    T Evaluate(const Intersection * si) const override;

    T Evaluate(vec2 uv) const override {

        ivec2 uvI(uv*vec2(float(_resU), float(_resV)));
        bool on = (uvI.x ^ uvI.y) & 1;
        return on ? _onColor : _offColor;
    }

    vec2 sample(TextureMapJacobian jacobian, const vec2 & uv,Float * pdf) const override {
        return vec2();
    }

    Float pdf(TextureMapJacobian jacobian, const vec2 & uv) const override {
        return 0;
    }

protected:
    T _onColor;
    T _offColor;
    int _resU;
    int _resV;
};

template < class T >
T CheckTexture < T >::Evaluate(const Intersection * si) const {
    return Evaluate(si->uv);
}


/**
 * @file 
 * @author JunPing Yuan
 * @brief 
 * @version 0.1
 * @date 2022/10/26
 *
 * @copyright Copyright (c) 2022
 *
 */
#include "Json.hpp"
#include "Common/Transform.hpp"

static vec3 randomOrtho(const vec3 &a)
{
    vec3 res;
    if (std::abs(a.x) > std::abs(a.y))
        res = vec3(0.0f, 1.0f, 0.0f);
    else
        res = vec3(1.0f, 0.0f, 0.0f);
    return normalize(cross(a,res));
}

static  mat4 rotYXZ(const vec3 &rot)
{
    vec3 r = rot*Constant::PI/180.0f;
    float c[] = {std::cos(r.x), std::cos(r.y), std::cos(r.z)};
    float s[] = {std::sin(r.x), std::sin(r.y), std::sin(r.z)};

    return mat4(
            c[1]*c[2] - s[1]*s[0]*s[2],   -c[1]*s[2] - s[1]*s[0]*c[2], -s[1]*c[0], 0.0f,
            c[0]*s[2],                     c[0]*c[2],      -s[0], 0.0f,
            s[1]*c[2] + c[1]*s[0]*s[2],   -s[1]*s[2] + c[1]*s[0]*c[2],  c[1]*c[0], 0.0f,
            0.0f,                          0.0f,       0.0f, 1.0f
    );
}

static void gramSchmidt(vec3 &a, vec3 &b, vec3 &c)
{
    a= normalize(a);
    b -= a*dot(a,b);
    if ( length2(b) < 1e-5)
        b = randomOrtho(a);
    else
        b= normalize(b);

    c -= a*dot(a,c);
    c -= b*dot(b,c);
    if ( length2(c) < 1e-5)
        c =cross(a,b);
    else
        c= normalize(c);
}

void glm::from_json(const Json & j,mat4 & transform) {
    if ( j.is_array() ) {   // given by explict matrix
        if ( j.size() != 16 ) {
            throw "Transfrom Matrix!=16";
        }
        for ( int i = 0 ; i < 16 ; i ++ )
            transform[i / 4][i % 4] = j[i];
    } else {

        vec3 x(1.0f, 0.0f, 0.0f);
        vec3 y(0.0f, 1.0f, 0.0f);
        vec3 z(0.0f, 0.0f, 1.0f);

        vec3 pos= getOptional(j,"position",vec3(0.0f));

        bool explicitX = false, explicitY = false, explicitZ = false;

        vec3 lookAt;
        if (containsAndGet(j,"look_at",lookAt) ) {
            z = lookAt - pos;
            explicitZ = true;
        }

        if ( containsAndGet(j,"up",y) ) {
            explicitY = true;
        }

//    explicitX = getField("x_axis", x) || explicitX;
//    explicitY = getField("y_axis", y) || explicitY;
//    explicitZ = getField("z_axis", z) || explicitZ;

        int id =
                ( explicitZ ? 4 : 0 ) +
                ( explicitY ? 2 : 0 ) +
                ( explicitX ? 1 : 0 );
        switch ( id ) {
            case 0:
                gramSchmidt(z, y, x);
                break;
            case 1:
                gramSchmidt(x, z, y);
                break;
            case 2:
                gramSchmidt(y, z, x);
                break;
            case 3:
                gramSchmidt(y, x, z);
                break;
            case 4:
                gramSchmidt(z, y, x);
                break;
            case 5:
                gramSchmidt(z, x, y);
                break;
            case 6:
                gramSchmidt(z, y, x);
                break;
            case 7:
                gramSchmidt(z, y, x);
                break;
        }

        if ( dot(cross(x, y), z) < 0.0f ) {
            if ( ! explicitX )
                x = -x;
            else if ( ! explicitY )
                y = - y;
            else
                z = - z;
        }

        vec3 scale;
        if ( j.contains("scale") ) {
            scale = j["scale"];
            x *= scale.x;
            y *= scale.y;
            z *= scale.z;
        }

        vec3 rot;
        if ( containsAndGet(j, "rotation", rot) ) {
            mat4 tform = rotYXZ(rot);
            x = mult(tform , vec4(x, 1.0));
            y = mult(tform , vec4(y, 1.0));
            z = mult(tform , vec4(z, 1.0));
        }
        transform = mat4(
                x[0], y[0], z[0], pos[0],
                x[1], y[1], z[1], pos[1],
                x[2], y[2], z[2], pos[2],
                0.0f, 0.0f, 0.0f, 1.0f
        );

    }  // end
}

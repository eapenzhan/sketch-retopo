#version 120
varying vec4 var_position;
varying vec3 var_normal;
varying vec4 var_color;
uniform vec3 uni_light_pos;
uniform int  uni_use_shading;

void main (void) {
    vec3 normal = normalize(var_normal);
    vec3 light  = normalize(uni_light_pos - var_position.xyz);
    gl_FragColor = var_color;
    if (uni_use_shading != 0) {
        float diffuse = max(dot(light, normal), 0);
        gl_FragColor.xyz *= diffuse;
    }
}

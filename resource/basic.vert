#version 120
varying vec4 var_position;
varying vec3 var_normal;
varying vec4 var_color;
uniform int  uni_symmetric_mirror;

void main(void) {
    vec4 vertex = gl_Vertex;
    vec3 normal = gl_Normal;
    
    if (uni_symmetric_mirror != 0) {
    	vertex.x *= -1.0;
    	normal.x *= -1.0;
    }
    
    var_position = gl_ModelViewMatrix * vertex;
    var_normal   = gl_NormalMatrix    * normal;
    var_color    = gl_Color;
    gl_Position  = gl_ModelViewProjectionMatrix * vertex;
}

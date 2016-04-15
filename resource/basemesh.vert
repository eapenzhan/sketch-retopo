#version 120
varying vec4 position;
varying vec3 normal;

void main(){
    gl_Position    = ftransform();
    gl_TexCoord[0] = gl_MultiTexCoord0;
    gl_TexCoord[1] = gl_MultiTexCoord1;
    position       = gl_Vertex;
    normal         = gl_NormalMatrix * gl_Normal;
}

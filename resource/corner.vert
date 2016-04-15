#version 120
uniform float fan_radius;
uniform int  uni_symmetric_mirror;

void main(){
    gl_FrontColor = gl_Color;
	vec3 fan_origin = gl_Vertex.xyz;
    vec3 fan_direction = gl_MultiTexCoord0.xyz;
    vec3 vertex = fan_origin + fan_radius * fan_direction + (fan_radius * 0.1) * gl_Normal;
    if (uni_symmetric_mirror != 0)
    	vertex.x *= -1.0;
    gl_Position = gl_ModelViewProjectionMatrix * vec4(vertex, 1.0);
}

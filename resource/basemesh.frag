#version 120
varying vec4 position;
varying vec3 normal;
uniform sampler2D material_texture;
uniform sampler2D texture_overlay;
uniform int symmetric;
uniform int line_mode;
uniform int show_expmap;
uniform float isoline_interval;
uniform float overlay_alpha;
uniform int overlay_v_flipped;

void main() {
    vec3 nnormal = normalize(normal);
    nnormal.y *= -1.0;
    vec2 material_texCoord = nnormal.xy * 0.45 + vec2(0.5, 0.5);
    vec4 fragColor = texture2D(material_texture, material_texCoord);
    if (symmetric != 0 && position[0] < 0.0)
        fragColor.rgb *= 0.7;
    
    if (show_expmap != 0) {
        vec2 uv = abs(gl_TexCoord[0].xy);
        if (uv != vec2(0.0)) {
            vec2 uv1 = uv / isoline_interval;
            vec2 uv2 = uv1 * 5.0;
            vec2 d1 = uv1 - floor(uv1);
            vec2 d2 = uv2 - floor(uv2);
            vec2 e1 = abs(d1 - vec2(0.5));
            vec2 e2 = abs(d2 - vec2(0.5));
            const float line_ratio = 0.03;
            if (e1.x > 0.5 - line_ratio || e1.y > 0.5 - line_ratio) fragColor.rgb = vec3(0.4);
            if (uv1.x < line_ratio) fragColor.rgb += vec3(-0.3,  0.5, -0.3);
            if (uv1.y < line_ratio) fragColor.rgb += vec3( 0.5, -0.3, -0.3);
            if (e2.x > 0.5 - line_ratio || e2.y > 0.5 - line_ratio) fragColor.rgb *= 0.75;
        }
    }
    vec2 uv = abs(gl_TexCoord[1].xy);
    if (overlay_v_flipped != 0)
        uv.y = 1.0 - uv.y;
    vec4 overlay_color = texture2D(texture_overlay, uv);
    float alpha = overlay_color.a * overlay_alpha;
    fragColor.rgb = (1 - alpha) * fragColor.rgb + alpha * overlay_color.rgb;
    
    if (line_mode != 0)
        fragColor.rgb *= 0.1;
    
    gl_FragColor = fragColor;
}

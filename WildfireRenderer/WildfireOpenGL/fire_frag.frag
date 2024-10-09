#version 330 core

// Shadertoy params
uniform float iTime;        // Time in seconds
uniform vec3 iResolution;   // Resolution of the viewport (width, height, depth)

out vec4 FragColor;

void main() {
    // Example usage of iTime and iResolution
    vec2 uv = gl_FragCoord.xy / iResolution.xy;  // Convert coordinates to 0..1 range
    vec3 color = 0.5 + 0.5 * cos(iTime + uv.xyx + vec3(0, 2, 4));
    
    FragColor = vec4(color, 1.0);  // Output final color
}
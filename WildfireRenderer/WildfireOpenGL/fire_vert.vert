#version 330 core
layout(location = 0) in vec2 aPos;       // Position attribute
layout(location = 1) in vec2 aTexCoord;  // Texture coordinate attribute

out vec2 TexCoord;  // Pass texture coordinate to the fragment shader

// Uniforms for iTime and iResolution (Optional if used here)
uniform float iTime;
uniform vec3 iResolution;

void main() {
    gl_Position = vec4(aPos, 0.0, 1.0); // Set the vertex position
    TexCoord = aTexCoord;               // Pass the texture coordinate
}
#version 330 core

layout(location = 0) in vec3 aPos;

// Uniforms for iTime and iResolution (Optional if used here)
uniform float iTime;
uniform vec3 iResolution;

void main() {
    gl_Position = vec4(aPos, 1.0);
}
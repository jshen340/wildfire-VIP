#version 330 core

uniform vec2 iResolution;
out vec2 fragCoord;

void main()
{
    // define a triangle in clip space
	vec2 vertices[3] = vec2[3](vec2(-1, -1), vec2(3, -1), vec2(-1, 3)); 
    
	// render area will be [-1, 1] x [-1, 1]
	gl_Position = vec4(vertices[gl_VertexID], 0, 1); 
    
	// convert to screen space [0, iResolution.x] x [0, iResolution.y]
	fragCoord = (0.5 * gl_Position.xy + vec2(0.5)) * iResolution;
}
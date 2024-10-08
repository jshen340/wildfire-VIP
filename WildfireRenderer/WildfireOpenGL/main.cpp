#include<iostream>
#include<glad/glad.h>
#include<GLFW/glfw3.h>

#pragma region Shader Source
/// <summary>
/// Vertex shader source, should not need to be edited since
/// shadertoy only exposes the fragment shader
/// </summary>
const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec2 aPos;

void main() {
    gl_Position = vec4(aPos, 0.0, 1.0);
}
)";

/// <summary>
/// Vertex shader source, should not need to be edited since
/// shadertoy only exposes the fragment shader
/// </summary>
const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;

uniform vec3 iResolution;
uniform float iTime;

void main() {
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    vec3 col = 0.5 + 0.5 * cos(iTime + uv.xyx + vec3(0, 2, 4));
    FragColor = vec4(col, 1.0);
}
)";
#pragma endregion Shader Source

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

unsigned int compileShader(unsigned int type, const char* source) {
	unsigned int shader = glCreateShader(type);
	glShaderSource(shader, 1, &source, NULL);
	glCompileShader(shader);

	int success;
	char infoLog[512];
	glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(shader, 512, NULL, infoLog);
		std::cerr << "ERROR::SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
	}

	return shader;
}

unsigned int createShaderProgram() {
	unsigned int vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
	unsigned int fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);

	unsigned int shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);

	int success;
	char infoLog[512];
	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
		std::cerr << "ERROR::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
	}

	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	return shaderProgram;
}

int main() {
	// Init GLFW
	if (!glfwInit()) {
		std::cerr << "Failed to initialize GLFW" << std::endl;
		return -1;
	}

    //Version 3.4 (10/8/2024)
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);

	//GLFW Core Profile (modern fucntions)
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	//Create new GLFWwindow popup, dimensions are currently 800 x 600 (same as shadertoy canvas)
	GLFWwindow* window = glfwCreateWindow(800, 600, "WildfireOpenGL", NULL, NULL);

	if (window == NULL) {
		std::cout << "GLFW window creation error" << std::endl;
		glfwTerminate();
		return -1;
	}
	//Link window to context
	glfwMakeContextCurrent(window);

	//Load OpenGL via GLAD
	gladLoadGL();

	//Viewport bounds in OpenGL & callback for window resizing
	glViewport(0, 0, 800, 600);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	//Compile & link shaders
	unsigned int shaderProgram = createShaderProgram();

	//Setup vertices before render loop -> renders quad that maps to screen dimensions
	float vertices[] = {
	-1.0f, -1.0f,  // Bottom-left
	 1.0f, -1.0f,  // Bottom-right
	 1.0f,  1.0f,  // Top-right
	-1.0f,  1.0f   // Top-left
	};

	unsigned int indices[] = {
		0, 1, 2,  // First triangle
		2, 3, 0   // Second triangle
	};

	unsigned int VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	// Bind VAO
	glBindVertexArray(VAO);

	// Bind VBO and load data
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// Bind EBO and load indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// Define vertex attribute (position)
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Unbind the VAO
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	//Main render loop
	while (!glfwWindowShouldClose(window)) {
		// Clear screen
		glClear(GL_COLOR_BUFFER_BIT);

		// Use the shader program
		glUseProgram(shaderProgram);

		// Set uniforms like time and resolution
		float timeValue = glfwGetTime();
		glUniform1f(glGetUniformLocation(shaderProgram, "iTime"), timeValue);

		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		glUniform3f(glGetUniformLocation(shaderProgram, "iResolution"), (float)width, (float)height, 1.0f);

		// Bind the VAO for the quad
		glBindVertexArray(VAO);

		// Draw the quad
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		// Swap buffers and poll for events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	//Terminate shader -> window -> GLFW -> program
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
	glDeleteBuffers(1, &EBO);
	glDeleteProgram(shaderProgram);
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
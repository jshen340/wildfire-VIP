#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <fstream>
#include <sstream>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define GL_VERTEX_SHADER 0x8B31
#define GL_FRAGMENT_SHADER 0x8B30

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

std::string readShaderFile(const std::string& filePath) {
    std::ifstream shaderFile;
    shaderFile.open(filePath);

    if (!shaderFile.is_open()) {
        std::cerr << "Failed to open shader file: " << filePath << std::endl;
        return "";
    }

    std::stringstream shaderStream;
    shaderStream << shaderFile.rdbuf();  // Read file's buffer contents into a string stream
    shaderFile.close();  // Close file handler

    return shaderStream.str();  // Convert stream into string
}

unsigned int compileShader(const std::string& vertexPath, const std::string& fragmentPath) {
    // Step 1: Read shader files
    std::string vertexCode = readShaderFile(vertexPath);
    std::string fragmentCode = readShaderFile(fragmentPath);

    if (vertexCode.empty() || fragmentCode.empty()) {
        std::cerr << "Shader code loading failed." << std::endl;
        return 0;  // Error handling
    }

    const char* vertexShaderSource = vertexCode.c_str();
    const char* fragmentShaderSource = fragmentCode.c_str();

    // Step 2: Compile Vertex Shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // Check compilation status
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    // Step 3: Compile Fragment Shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // Check compilation status
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    // Step 4: Link Shaders
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // Check linking status
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }

    // Step 5: Clean up
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}
void errorCallback(int error, const char* description) {
    std::cerr << "GLFW Error: " << description << std::endl;
}

int main() {
    // Init GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Version 3.4 (10/8/2024)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    // GLFW Core Profile (modern functions)
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);  // Ensure the window will be visible


    // Create new GLFWwindow popup, dimensions are currently 800 x 600
    GLFWwindow* window = glfwCreateWindow(800, 600, "Texture Rendering OpenGL", NULL, NULL);

    if (window == NULL) {
        std::cout << "GLFW window creation error" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Link window to context
    glfwMakeContextCurrent(window);

    // Load OpenGL via GLAD
    gladLoadGL();

    // Viewport bounds in OpenGL & callback for window resizing
    glViewport(0, 0, 800, 600);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Compile & link shaders
    unsigned int shaderProgram = compileShader("fire_vert.vert", "fire_frag.frag");

    // Load texture
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    // Set texture wrapping and filtering options
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Load image
    int texWidth, texHeight, texChannels;
    unsigned char* data = stbi_load("IMG_1159.jpg", &texWidth, &texHeight, &texChannels, 0);
    if (data) {
        // Generate the texture
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texWidth, texHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    else {
        std::cerr << "Failed to load texture" << std::endl;
    }
    stbi_image_free(data);

    // Setup vertices and texture coordinates
    float vertices[] = {
        // positions        // texture coordinates
        -1.0f, -1.0f,      0.0f, 0.0f,  // Bottom-left
         1.0f, -1.0f,      1.0f, 0.0f,  // Bottom-right
         1.0f,  1.0f,      1.0f, 1.0f,  // Top-right
        -1.0f,  1.0f,      0.0f, 1.0f   // Top-left
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
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Define texture coordinate attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Unbind the VAO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Main render loop
    while (!glfwWindowShouldClose(window)) {
        // Clear screen
        glClear(GL_COLOR_BUFFER_BIT);

        // Use the shader program
        glUseProgram(shaderProgram);

        // Bind the texture
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glUniform1i(glGetUniformLocation(shaderProgram, "ourTexture"), 0);

        // Bind the VAO for the quad
        glBindVertexArray(VAO);

        // Draw the quad
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // Swap buffers and poll for events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Clean up
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

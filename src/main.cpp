#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <random>
#include <stb_image_write.h>

#include "Common.h"
#include "MPM.h"

const int d = 2;
const real dt = 1e-4;

using VectorD = Vector<real, d>;
using VectorDi = Vector<int, d>;

void framebuffer_size_callback(GLFWwindow *window, int width, int height);
void processInput(GLFWwindow *window);
void Three_Block_Init();
void saveImage(char *filepath, GLFWwindow *w);
MPM<d> mpm;

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 800;

const char *vertexShaderSource = "#version 410 core\n"
                                 "layout (location = 0) in vec3 aPos;\n"
                                 "void main()\n"
                                 "{\n"
                                 "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
                                 "   gl_PointSize = 5.0;\n"
                                 "}\0";

const char *fragmentShaderSource = "#version 410 core\n"
                                   "out vec4 FragColor;\n"
                                   "void main()\n"
                                   "{\n"
                                   "   FragColor = vec4(1);\n"
                                   "}\n\0";

const char *vertexPlaneShaderSource = "#version 410 core\n"
                                      "layout (location = 0) in vec3 aPos;\n"
                                      "void main()\n"
                                      "{\n"
                                      "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
                                      "}\0";

const char *fragmentPlaneShaderSource = "#version 410 core\n"
                                        "out vec4 FragColorPlane;\n"
                                        "void main()\n"
                                        "{\n"
                                        "   FragColorPlane = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
                                        "}\n\0";

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "MPM", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    unsigned int vertexPlaneShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexPlaneShader, 1, &vertexPlaneShaderSource, NULL);
    glCompileShader(vertexPlaneShader);
    // check for shader compile errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexPlaneShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexPlaneShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
    }
    // fragment shader
    unsigned int fragmentPlaneShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentPlaneShader, 1, &fragmentPlaneShaderSource, NULL);
    glCompileShader(fragmentPlaneShader);
    // check for shader compile errors
    glGetShaderiv(fragmentPlaneShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentPlaneShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
    }
    // link shaders
    unsigned int shaderProgramPlane = glCreateProgram();
    glAttachShader(shaderProgramPlane, vertexPlaneShader);
    glAttachShader(shaderProgramPlane, fragmentPlaneShader);
    glLinkProgram(shaderProgramPlane);
    // check for linking errors
    glGetProgramiv(shaderProgramPlane, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(shaderProgramPlane, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n"
                  << infoLog << std::endl;
    }
    glDeleteShader(vertexPlaneShader);
    glDeleteShader(fragmentPlaneShader);

    float vertices[] = {
        0.99f, 0.99f,
        0.99f, 0.0f,
        0.99f, 0.0f,
        0.0f, 0.0f,
        0.0f, 0.0f,
        0.0f, 0.99f,
        0.0f, 0.99f,
        0.99f, 0.99f};

    unsigned int VBO_plane, VAO_plane;
    glGenVertexArrays(1, &VAO_plane);
    glGenBuffers(1, &VBO_plane);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO_plane);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_plane);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, d, GL_FLOAT, GL_FALSE, d * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
    glBindBuffer(GL_ARRAY_BUFFER, VBO_plane);

    // build and compile our shader program
    // ------------------------------------
    // vertex shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
    }
    // fragment shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    // check for shader compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
    }
    // link shaders
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    // check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n"
                  << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    // add a new set of vertices to form a second triangle (a total of 6 vertices); the vertex attribute configuration remains the same (still one 3-float position vector per vertex)
    Three_Block_Init();
    mpm.Initialize();

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, mpm.particles.Size() * sizeof(mpm.particles.x->data()), mpm.particles.x->data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, d, GL_FLOAT, GL_FALSE, d * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    // render loop
    // -----------
    int step = 0;
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgramPlane);
        glBindVertexArray(VAO_plane); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
        glDrawArrays(GL_LINES, 0, 8);

        glUseProgram(shaderProgram);
        for (int i = 0; i < 5; i++)
        {
            mpm.Advance(dt);
        }
        char fname[128];
        sprintf(fname, "./image/test/frame_%05d.png", step);
        saveImage(fname, window);
        glBindVertexArray(VAO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, mpm.particles.Size() * sizeof(mpm.particles.x->data()), mpm.particles.x->data());
        glDrawArrays(GL_POINTS, 0, mpm.particles.Size());
        step++;
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

void Add_Particle(VectorD &pos, VectorD &u)
{
    int i = mpm.particles.Add_Element(); ////return the last element's index
    mpm.particles.X(i) = pos;
    mpm.particles.M(i) = 1.;
    mpm.particles.V(i) = u;
}

void Add_Object(VectorD &center, int num_samplings, VectorD &u)
{
    for (int i = 0; i < num_samplings; i++)
    {
        VectorD pos = 0.08 * (VectorD::Random() - VectorD::Ones()) + center;
        Add_Particle(pos, u);
    }
}

void Three_Block_Init()
{
    real dx = .015;
    VectorDi points_range = 15 * VectorDi::Ones();
    VectorD u = VectorD::Zero();

    VectorD center1 = VectorD(0.4, 0.6) + points_range.template cast<real>() * dx / 2;
    VectorD center2 = VectorD(0.52, 0.8) + points_range.template cast<real>() * dx / 2;
    VectorD center3 = VectorD(0.54, 0.4) + points_range.template cast<real>() * dx / 2;

    Add_Object(center1, points_range.prod(), u);
    Add_Object(center2, points_range.prod(), u);
    Add_Object(center3, points_range.prod(), u);
}

void saveImage(char *filepath, GLFWwindow *w)
{
    int width, height;
    glfwGetFramebufferSize(w, &width, &height);
    GLsizei nrChannels = 3;
    GLsizei stride = nrChannels * width;
    stride += (stride % 4) ? (4 - stride % 4) : 0;
    GLsizei bufferSize = stride * height;
    std::vector<char> buffer(bufferSize);
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(400, 400, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
    stbi_flip_vertically_on_write(true);
    stbi_write_png(filepath, width, height, nrChannels, buffer.data(), stride);
}

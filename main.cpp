#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
//
#include <glbinding/gl/gl.h>
#include <glbinding/glbinding.h>
using namespace gl;
//
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
//
#include <glm/glm.hpp>
//
#include <glm/ext.hpp>

float width = 800;
float height = 450;
glm::vec2 fov{M_PI_4 * width / height, M_PI_4};
glm::vec3 origin{};
glm::vec3 up{0, 1, 0};
glm::vec3 camera{10, 0.0f, M_PI_2};

int main(void) {
  using namespace std;

  // Generate random points in sphere.
  mt19937 rng{random_device{}()};
  uniform_real_distribution<float> dist{-1.0f, 1.0f};

  const size_t samples = 1000;
  vector<glm::vec3> points(samples);
  for (auto& p : points) p = glm::vec3{dist(rng), dist(rng), dist(rng)};

  glfwSetErrorCallback([](int error, const char* description) {
    throw runtime_error{"GLFW Error " + to_string(error) + ": " + description};
  });

  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_SAMPLES, 4);

  auto window =
      glfwCreateWindow(width, height, "Ray Marcher", nullptr, nullptr);
  glfwMakeContextCurrent(window);
  glbinding::initialize(glfwGetProcAddress);

  glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int w, int h) {
    width = w;
    height = h;
    fov.x = fov.y * width / height;
  });
  glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode,
                                int action, int mods) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
      glfwSetWindowShouldClose(window, GLFW_TRUE);
  });
  glfwSetScrollCallback(window, [](GLFWwindow* window, double x, double y) {
    camera.x *= exp(-0.1f * float(y));
    // camera = exp(-0.1f * float(y)) * (camera - origin) + origin;
  });

  GLuint vertex_array;
  glGenVertexArrays(1, &vertex_array);
  glBindVertexArray(vertex_array);

  // GLuint element_buffer;
  // glGenBuffers(1, &element_buffer);
  // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer);
  // glBufferData(GL_ELEMENT_ARRAY_BUFFER, elements.size() * sizeof(uint32_t),
  //              elements.data(), GL_STATIC_DRAW);

  GLuint vertex_buffer;
  glGenBuffers(1, &vertex_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * points.size(),
               points.data(), GL_STATIC_DRAW);

  auto vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  const char* vertex_shader_text =
      "#version 330\n"
      "uniform mat4 MVP;"
      "attribute vec3 vPos;"
      "void main() {"
      "  gl_Position = MVP * vec4(vPos, 1.0);"
      "}";
  glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
  glCompileShader(vertex_shader);
  auto fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  const char* fragment_shader_text =
      "#version 330\n"
      "void main() {"
      "  gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);"
      "}";
  glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
  glCompileShader(fragment_shader);
  auto program = glCreateProgram();
  glAttachShader(program, vertex_shader);
  glAttachShader(program, fragment_shader);
  glLinkProgram(program);
  auto mvp_location = glGetUniformLocation(program, "MVP");
  auto vpos_location = glGetAttribLocation(program, "vPos");

  glEnableVertexAttribArray(vpos_location);
  glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3),
                        (void*)0);

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glPointSize(5.0f);
  glEnable(GL_MULTISAMPLE);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_POINT_SPRITE);
  glEnable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glm::vec2 old_mouse_pos{};
  glm::vec2 mouse_pos{};

  while (!glfwWindowShouldClose(window)) {
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    old_mouse_pos = mouse_pos;
    mouse_pos = glm::vec2(xpos, ypos);

    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS) {
      const auto mouse_move = mouse_pos - old_mouse_pos;
      camera.z += mouse_move.x * 0.01;
      camera.y += mouse_move.y * 0.01;
      const constexpr float eye_altitude_max_abs = M_PI_2 - 0.0001f;
      camera.y = clamp(camera.y, -eye_altitude_max_abs, eye_altitude_max_abs);
    }

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glm::mat4x4 m{1.0f};
    // m = rotate(m, (float)glfwGetTime(), glm::vec3(1, 1, 1));
    const auto v = glm::lookAt(
        origin + camera.x * glm::vec3{cos(camera.y) * cos(camera.z),  //
                                      sin(camera.y),                  //
                                      cos(camera.y) * sin(camera.z)},
        origin, up);
    const auto p = glm::perspective(fov.y, width / height, 0.1f, 100.f);
    glm::mat4 mvp = p * v * m;

    glUseProgram(program);
    glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));

    // glPointSize(1.0f);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * points.size(),
                 points.data(), GL_DYNAMIC_DRAW);
    glDrawArrays(GL_POINTS, 0, points.size());
    // glDrawElements(GL_TRIANGLES, elements.size(), GL_UNSIGNED_INT, nullptr);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwDestroyWindow(window);
  glfwTerminate();
}
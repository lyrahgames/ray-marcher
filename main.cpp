#include <cmath>
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
#include <glm/gtx/norm.hpp>

float width = 800;
float height = 450;
glm::vec2 fov{M_PI_4 * width / height, M_PI_4};
glm::vec3 origin{};
glm::vec3 up{0, 1, 0};
glm::vec3 camera{5, 0.0f, M_PI_2};

std::vector<unsigned char> data{};

void update() {
  // Camera vectors
  auto direction =
      -camera.x * glm::vec3{cos(camera.y) * cos(camera.z), sin(camera.y),
                            cos(camera.y) * sin(camera.z)};
  auto position = origin - direction;
  direction = normalize(direction);
  auto right = cross(direction, up);
  right = normalize(right);
  auto cup = cross(right, direction);
  float pixel_size = 2.0f * tan(0.5f * fov.y) / height;

  data.resize(3 * width * height);

#pragma omp parallel for
  for (size_t j = 0; j < height; ++j) {
    for (size_t i = 0; i < width; ++i) {
      // Construct primary ray.
      auto ray_o = position;
      auto ray_d = direction + (pixel_size * (i - width / 2)) * right +
                   (pixel_size * (j - height / 2)) * cup;
      ray_d = normalize(ray_d);

      float total = 0.0f;
      int steps = 0;
      float last = INFINITY;
      const int max_steps = 100;
      for (; steps < max_steps; ++steps) {
        auto p = ray_o + total * ray_d;
        //
        // float distance = std::abs(glm::length(p - origin) - 1.0f);
        auto z = p;
        float dr = 1.0;
        float r = 0.0;
        for (int i = 0; i < 1000; i++) {
          r = length(z);
          if (r > 4) break;
          // convert to polar coordinates
          float theta = std::acos(z.z / r);
          float phi = std::atan(z.y / z.x);
          const float Power = 8;
          dr = std::pow(r, Power - 1.0) * Power * dr + 1.0;
          // scale and rotate the point
          float zr = pow(r, Power);
          theta = theta * Power;
          phi = phi * Power;
          // convert back to cartesian coordinates
          z = zr * glm::vec3(std::sin(theta) * std::cos(phi),
                             std::sin(phi) * std::sin(theta), std::cos(theta));
          z += p;
        }
        float distance = 0.5 * std::log(r) * r / dr;
        //
        total += distance;
        if (distance < 1e-3f) break;
        last = distance;
      }
      float color = (1.0f - float(steps) / float(max_steps));
      // float color = 0.5f * (dot(ray_d, direction) + 1.0f);

      const auto index = 3 * (j * width + i);
      data[index + 0] = 255.0f * color;
      data[index + 1] = 255.0f * color;
      data[index + 2] = 255.0f * color;
    }
  }
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB,
               GL_UNSIGNED_BYTE, data.data());
}

int main(void) {
  using namespace std;

  // Generate random points in sphere.
  mt19937 rng{random_device{}()};
  uniform_real_distribution<float> dist{-1.0f, 1.0f};

  vector<glm::vec2> points{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  vector<uint32_t> elements{0, 1, 2, 0, 2, 3};

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

  GLuint element_buffer;
  glGenBuffers(1, &element_buffer);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, elements.size() * sizeof(uint32_t),
               elements.data(), GL_STATIC_DRAW);

  GLuint vertex_buffer;
  glGenBuffers(1, &vertex_buffer);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
  glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2) * points.size(),
               points.data(), GL_STATIC_DRAW);

  auto vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  const char* vertex_shader_text =
      "#version 330\n"
      "uniform mat4 MVP;"
      "attribute vec2 vPos;"
      "out vec2 texuv;"
      "void main() {"
      "  gl_Position = MVP * vec4(vPos, 0.0, 1.0);"
      "  texuv = vPos;"
      "}";
  glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
  glCompileShader(vertex_shader);
  auto fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  const char* fragment_shader_text =
      "#version 330\n"
      "in vec2 texuv;"
      "uniform sampler2D tex;"
      "void main() {"
      "  gl_FragColor = texture(tex, texuv);"
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
  glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2),
                        (void*)0);

  GLuint texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB,
  //              GL_UNSIGNED_BYTE, data.data());
  // glGenerateMipmap(GL_TEXTURE_2D);

  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
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

    update();

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // glm::mat4x4 m{1.0f};
    // // m = rotate(m, (float)glfwGetTime(), glm::vec3(1, 1, 1));
    // const auto v = glm::lookAt(
    //     origin + camera.x * glm::vec3{cos(camera.y) * cos(camera.z),  //
    //                                   sin(camera.y),                  //
    //                                   cos(camera.y) * sin(camera.z)},
    //     origin, up);
    // const auto p = glm::perspective(fov.y, width / height, 0.1f, 100.f);
    // glm::mat4 mvp = p * v * m;

    glm::mat4 mvp = glm::ortho(0, 1, 0, 1, 0, 1);

    glUseProgram(program);
    glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(mvp));

    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * points.size(),
                 points.data(), GL_DYNAMIC_DRAW);
    // glDrawArrays(GL_POINTS, 0, points.size());
    glDrawElements(GL_TRIANGLES, elements.size(), GL_UNSIGNED_INT, nullptr);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwDestroyWindow(window);
  glfwTerminate();
}
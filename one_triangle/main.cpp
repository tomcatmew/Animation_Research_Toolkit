// Modified by Yifei Chen for OpenGL review and learning

#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cstdio>

static void error_callback(int error, const char* description)
{
  fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);
}

int main(void)
{
  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(840, 680, "Assignment_1_One_Triangle", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);

  while (!glfwWindowShouldClose(window))
  {
    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
    glBegin(GL_TRIANGLES);
    glColor3f(0.7f, 1.f, 0.f);
    glVertex3f(-0.6f, -0.4f, 0.f);

    glColor3f(0.7f, 1.f, 0.f);
    glVertex3f(0.6f, -0.4f, 0.f);

    glColor3f(0.1f, 1.f, 1.f);
    glVertex3f(0.0f, 0.8f, 0.f);

    glEnd();

    // give triangle a box !
    glBegin(GL_LINE_STRIP);
    glColor3f(1.f, 1.f, 1.f);
    glVertex3f(-0.9f, -0.8f, 0.f);
    glVertex3f(0.9f, -0.8f, 0.f);
    glVertex3f(0.9f, 0.9f, 0.f);
    glVertex3f(-0.9f, 0.9f, 0.f);
    glVertex3f(-0.9f, -0.8f, 0.f);
    glEnd();

    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
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
  double origin[] = { -0.9f, -0.5f, 0.f };
  double line_1[] = { 0.9f, -0.5f, 0.f };
  double line_2[] = { -0.9f, 0.5f, 0.f };
  double t1[] = { line_1[0] - origin[0],line_1[1] - origin[1] ,line_1[2] - origin[2] };
  double t2[] = { line_2[0] - origin[0],line_2[1] - origin[1] ,line_2[2] - origin[2] };

  double line_3[] = { t1[0] + t2[0] + origin[0],t1[1] + t2[1] + origin[1],t1[2] + t2[2] + origin[2] };

  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
    exit(EXIT_FAILURE);
  window = glfwCreateWindow(1340, 880, "Assignment_1_One_Triangle", NULL, NULL);
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
    double scale_ratio = 6.0f;
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-ratio * scale_ratio, ratio * scale_ratio, -1.f * scale_ratio, 1.f * scale_ratio, 100.f * scale_ratio, -100.f * scale_ratio);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
    //glBegin(GL_TRIANGLES);
    //glColor3f(0.7f, 1.f, 0.f);
    //glVertex3f(-0.6f, -0.4f, 0.f);

    //glColor3f(0.7f, 1.f, 0.f);
    //glVertex3f(0.6f, -0.4f, 0.f);

    //glColor3f(0.1f, 1.f, 1.f);
    //glVertex3f(0.0f, 0.8f, 0.f);

    //glEnd();

    // give triangle a box !
    glBegin(GL_LINE_STRIP);
    glColor3f(1.f, 1.f, 1.f);
    glVertex3f(origin[0], origin[1], origin[2]);
    glVertex3f(line_1[0], line_1[1], line_1[2]);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3f(1.f, 1.f, 1.f);
    glVertex3f(origin[0], origin[1], origin[2]);
    glVertex3f(line_2[0], line_2[1], line_2[2]);
    //glVertex3f(0.9f, 0.9f, 0.f);
    //glVertex3f(-0.9f, -0.8f, 0.f);
    glEnd();
    glBegin(GL_LINE_STRIP);
    glColor3f(1.f, 1.f, 1.f);
    glVertex3f(origin[0], origin[1], origin[2]);
    glVertex3f(line_3[0], line_3[1], line_3[2]);
    glEnd();

    glfwSwapBuffers(window);
    glfwPollEvents();
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
// DH2323 skeleton code, Lab3 (SDL2 version)
#include "SDL2Auxiliary.h"
#include "TestModel.h"
#include <algorithm> //for max()
#include <cstddef>
#include <cstdio>
#include <glm/glm.hpp>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

using namespace std;
using glm::mat3;
using glm::vec3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL2Aux *sdlAux;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001);
mat3 R(1.0f);
float yaw = 0;
const float moveSpeed = 0.05f;
vec3 currentColor(0, 0, 0);

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update(void);
void Draw(void);

void VertexShader(const vec3 &v, glm::ivec2 &p);
void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2> &result);
void DrawLineSDL(glm::ivec2 a, glm::ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3> &vertices);
void ComputePolygonRows(const vector<glm::ivec2> &vertexPixels,
                        vector<glm::ivec2> &leftPixels,
                        vector<glm::ivec2> &rightPixels);
void DrawPolygonRows(const vector<glm::ivec2> &leftPixels,
                     const vector<glm::ivec2> &rightPixels);
void DrawPolygon(const vector<vec3> &vertices);

int main(int argc, char *argv[]) {
  LoadTestModel(triangles); // Load model
  sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
  t = SDL_GetTicks(); // Set start value for timer.

  while (!sdlAux->quitEvent()) {
    Update();
    Draw();
  }

  sdlAux->saveBMP("screenshot.bmp");
  return 0;
}

void Update(void) {
  // Compute frame time:
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t);
  t = t2;
  cout << "Render time: " << dt << " ms." << endl;

  // int dx, dy;
  // SDL_GetRelativeMouseState(&dx, &dy);
  //
  // float mouseSensitivity = 0.002f;
  // yaw -= dx * mouseSensitivity;

  const Uint8 *keystate = SDL_GetKeyboardState(NULL);
  if (keystate[SDL_SCANCODE_UP]) {
    // Move camera forward
    cameraPos += moveSpeed * vec3(R[0][2], R[1][2], R[2][2]);
  }
  if (keystate[SDL_SCANCODE_DOWN]) {
    // Move camera backward
    cameraPos -= moveSpeed * vec3(R[0][2], R[1][2], R[2][2]);
  }
  if (keystate[SDL_SCANCODE_LEFT]) {
    // Move camera to the left
    yaw -= moveSpeed;
  }
  if (keystate[SDL_SCANCODE_RIGHT]) {
    // Move camera to the right
    yaw += moveSpeed;
  }
  if (keystate[SDL_SCANCODE_W]) {
  }
  if (keystate[SDL_SCANCODE_S]) {
  }
  if (keystate[SDL_SCANCODE_A]) {
  }
  if (keystate[SDL_SCANCODE_D]) {
  }
  if (keystate[SDL_SCANCODE_Q]) {
  }
  if (keystate[SDL_SCANCODE_E]) {
  }
  R = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
}

void Draw() {
  sdlAux->clearPixels();

  for (int i = 0; i < triangles.size(); ++i) {
    vector<vec3> vertices(3);

    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    currentColor = triangles[i].color;

    // DrawPolygonEdges(vertices);
    DrawPolygon(vertices);
  }

  sdlAux->render();
}

void VertexShader(const vec3 &v, glm::ivec2 &p) {
  float f = SCREEN_HEIGHT;
  vec3 v1 = R * (v - cameraPos);
  p.x = f * (v1.x / v1.z) + (SCREEN_WIDTH / 2.f);
  p.y = f * (v1.y / v1.z) + (SCREEN_HEIGHT / 2.f);
}

void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2> &result) {
  int N = result.size();
  glm::vec2 step = glm::vec2(b - a) / float(glm::max(N - 1, 1));
  glm::vec2 current(a);
  for (int i = 0; i < N; ++i) {
    result[i] = current;
    current += step;
  }
}

void DrawLineSDL(glm::ivec2 a, glm::ivec2 b, vec3 color) {
  glm::ivec2 delta = glm::abs(a - b);
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<glm::ivec2> line(pixels);
  Interpolate(a, b, line);
  for (glm::ivec2 pixel : line) {
    sdlAux->putPixel(pixel.x, pixel.y, color);
  }
}

void DrawPolygonEdges(const vector<vec3> &vertices) {
  int V = vertices.size();

  // Transform each vertex from 3D world position to 2D image position
  vector<glm::ivec2> projectedVertices(V);
  for (int i = 0; i < V; ++i) {
    VertexShader(vertices[i], projectedVertices[i]);
  }

  // Loop over all vertices and draw the edge from it to the next vertex
  for (int i = 0; i < V; ++i) {
    int j = (i + 1) % V;
    vec3 color(1, 1, 1);
    DrawLineSDL(projectedVertices[i], projectedVertices[j], color);
  }
}

void ComputePolygonRows(const vector<glm::ivec2> &vertexPixels,
                        vector<glm::ivec2> &leftPixels,
                        vector<glm::ivec2> &rightPixels) {
  auto minmaxY = std::minmax_element(
      vertexPixels.begin(), vertexPixels.end(),
      [](const glm::ivec2 &a, const glm::ivec2 &b) { return a.y < b.y; });

  auto minY = minmaxY.first->y;
  auto maxY = minmaxY.second->y;
  // printf("minY: %d, maxY: %d\n", minY, maxY);

  auto rows = maxY - minY + 1;

  leftPixels.resize(rows);
  rightPixels.resize(rows);
  for (size_t i = 0; i < rows; ++i) {
    leftPixels[i].x = numeric_limits<int>::max();
    leftPixels[i].y = minY + i;
    rightPixels[i].x = numeric_limits<int>::min();
    rightPixels[i].y = minY + i;
  }

  for (int i = 0; i < vertexPixels.size(); ++i) {

    glm::ivec2 a = vertexPixels[i];
    glm::ivec2 b = vertexPixels[(i + 1) % vertexPixels.size()];

    int deltaY = glm::abs(b.y - a.y);
    int pixels = deltaY + 1;

    vector<glm::ivec2> edge(pixels);
    Interpolate(a, b, edge);

    for (const glm::ivec2 &p : edge) {
      int y = p.y - minY;
      if (y < leftPixels.size()) {
        if (p.x < leftPixels[y].x)
          leftPixels[y].x = p.x;
        if (p.x > rightPixels[y].x)
          rightPixels[y].x = p.x;
      }
    }
  }
}

void DrawPolygonRows(const vector<glm::ivec2> &leftPixels,
                     const vector<glm::ivec2> &rightPixels) {
  for (size_t i = 0; i < leftPixels.size(); ++i) {
    auto y = leftPixels[i].y;
    auto x_start = rightPixels[i].x;
    auto x_end = leftPixels[i].x;
    if (x_start > x_end)
      std::swap(x_start, x_end);
    for (auto x = x_start; x < x_end; ++x) {
      sdlAux->putPixel(x, y, currentColor);
    }
  }
}

void DrawPolygon(const vector<vec3> &vertices) {
  int V = vertices.size();
  vector<glm::ivec2> vertexPixels(V);
  for (int i = 0; i < V; ++i) {
    VertexShader(vertices[i], vertexPixels[i]);
  }
  vector<glm::ivec2> leftPixels;
  vector<glm::ivec2> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(leftPixels, rightPixels);
}

// DH2323 skeleton code, Lab3 (SDL2 version)
#include "SDL2Auxiliary.h"
#include "TestModel.h"
#include <algorithm> //for max()
#include <glm/glm.hpp>
#include <iostream>

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

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update(void);
void Draw(void);

void VertexShader(const vec3 &v, glm::ivec2 &p);
void Interpolate(glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2> &result);
void DrawLineSDL(glm::ivec2 a, glm::ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3> &vertices);

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

  const Uint8 *keystate = SDL_GetKeyboardState(NULL);
  if (keystate[SDL_SCANCODE_UP]) {
    // Move camera forward
  }
  if (keystate[SDL_SCANCODE_DOWN]) {
    // Move camera backward
  }
  if (keystate[SDL_SCANCODE_LEFT]) {
    // Move camera to the left
  }
  if (keystate[SDL_SCANCODE_RIGHT]) {
    // Move camera to the right
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
}

void Draw() {
  sdlAux->clearPixels();

  for (int i = 0; i < triangles.size(); ++i) {
    vector<vec3> vertices(3);

    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    // Add drawing
    for (int v = 0; v < 3; ++v) {
      glm::ivec2 projPos;
      vertices[v] = vertices[v] - cameraPos;
      VertexShader(vertices[v], projPos);
      vec3 color(1, 1, 1);
      sdlAux->putPixel(projPos.x, projPos.y, color);
    }
  }

  sdlAux->render();
}

void VertexShader(const vec3 &v, glm::ivec2 &p) {
  float f = SCREEN_HEIGHT;
  p.x = f * (v.x / v.z) + (SCREEN_WIDTH / 2.f);
  p.y = f * (v.y / v.z) + (SCREEN_HEIGHT / 2.f);
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

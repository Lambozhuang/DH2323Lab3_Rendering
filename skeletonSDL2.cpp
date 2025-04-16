// DH2323 skeleton code, Lab3 (SDL2 version)
#include "SDL2Auxiliary.h"
#include "TestModel.h"
#include <algorithm> //for max()
#include <cstddef>
#include <cstdio>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

using namespace std;
using glm::mat3;
using glm::vec2;
using glm::vec3;

struct Vertex {
  vec3 position;
};

struct Pixel {
  int x, y;
  float zinv;
  vec3 pos3d;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL2Aux *sdlAux;
int t;
vector<Triangle> triangles;

vec3 currentColor(0, 0, 0);
vec3 currentNormal;
vec3 currentReflectance;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

// Camera model
vec3 cameraPos(0, 0, -3.001);
mat3 R(1.0f);
float yaw = 0;
const float moveSpeed = 0.05f;

// Light model
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 11.f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1, 1, 1);

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update(void);
void Draw(void);

void VertexShader(const Vertex &v, Pixel &p);
void PixelShader(const Pixel &p);
void Interpolate(Pixel a, Pixel b, vector<Pixel> &result);
void DrawLineSDL(Pixel a, Pixel b, vec3 color);
void DrawPolygonEdges(const vector<Vertex> &vertices);
void ComputePolygonRows(const vector<Pixel> &vertexPixels,
                        vector<Pixel> &leftPixels, vector<Pixel> &rightPixels);
void DrawPolygonRows(const vector<Pixel> &leftPixels,
                     const vector<Pixel> &rightPixels);
void DrawPolygon(const vector<Vertex> &vertices);

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

  // Clear depth buffer
  for (int i = 0; i < SCREEN_HEIGHT; ++i) {
    for (int j = 0; j < SCREEN_WIDTH; ++j) {
      depthBuffer[i][j] = 0;
    }
  }

  for (int i = 0; i < triangles.size(); ++i) {
    vector<Vertex> vertices(3);

    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;

    currentNormal = triangles[i].normal;
    currentReflectance = vec3(1, 1, 1);
    currentColor = triangles[i].color;

    // DrawPolygonEdges(vertices);
    DrawPolygon(vertices);
  }

  sdlAux->render();
}

// Given a 3D vertex, transform it to 2D screen coordinates
void VertexShader(const Vertex &v, Pixel &p) {
  float f = SCREEN_HEIGHT;
  vec3 v1 = R * (v.position - cameraPos);
  p.x = f * (v1.x / v1.z) + (SCREEN_WIDTH / 2.f);
  p.y = f * (v1.y / v1.z) + (SCREEN_HEIGHT / 2.f);
  p.zinv = 1.f / v1.z;
  p.pos3d = v.position;
}

void PixelShader(const Pixel &p) {
  int x = p.x;
  int y = p.y;
  if (x < 0 || x >= SCREEN_WIDTH || y < 0 || y >= SCREEN_HEIGHT) {
    return;
  }
  // Check if the pixel is closer than the current depth value
  if (depthBuffer[y][x] < p.zinv) {
    depthBuffer[y][x] = p.zinv;
    // Compute the color based on the light model
    vec3 lightDirection = glm::normalize(lightPos - p.pos3d);
    float distanceToLight = glm::distance(lightPos, p.pos3d);
    vec3 illumination = (lightPower * glm::max(0.f, glm::dot(lightDirection, currentNormal))) / (4 * distanceToLight * distanceToLight * glm::pi<float>()) + indirectLightPowerPerArea;
    // Draw the pixel with the current color
    sdlAux->putPixel(x, y, illumination * currentColor);
  }
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result) {
  int N = result.size();

  vec2 stepXY = vec2(b.x - a.x, b.y - a.y) / float(glm::max(N - 1, 1));
  float stepZinv = (b.zinv - a.zinv) / float(glm::max(N - 1, 1));
  vec3 stepPos3d = (b.pos3d - a.pos3d) / float(glm::max(N - 1, 1));

  float currentX = a.x;
  float currentY = a.y;
  float currentZinv = a.zinv;
  vec3 currentPos3d = a.pos3d;

  for (int i = 0; i < N; ++i) {
    
    // Should be careful with float -> int
    result[i].x = static_cast<int>(round(currentX)); // Round instead of truncate
    result[i].y = static_cast<int>(round(currentY)); // Round instead of truncate
    result[i].zinv = currentZinv;
    result[i].pos3d = currentPos3d;

    currentX += stepXY.x;
    currentY += stepXY.y;
    currentZinv += stepZinv;
    currentPos3d += stepPos3d;
  }
}

void DrawLineSDL(Pixel a, Pixel b, vec3 color) {
  glm::ivec2 delta = glm::ivec2(glm::abs(a.x - b.x), glm::abs(a.y - b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);
  for (Pixel pixel : line) {
    sdlAux->putPixel(pixel.x, pixel.y, color);
  }
}

void DrawPolygonEdges(const vector<Vertex> &vertices) {
  int V = vertices.size();

  // Transform each vertex from 3D world position to 2D image position
  vector<Pixel> projectedVertices(V);
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

void ComputePolygonRows(const vector<Pixel> &vertexPixels,
                        vector<Pixel> &leftPixels, vector<Pixel> &rightPixels) {
  auto minmaxY = std::minmax_element(
      vertexPixels.begin(), vertexPixels.end(),
      [](const Pixel &a, const Pixel &b) { return a.y < b.y; });

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

    Pixel a = vertexPixels[i];
    Pixel b = vertexPixels[(i + 1) % vertexPixels.size()];

    int deltaY = glm::abs(b.y - a.y);
    int pixels = deltaY + 1;

    vector<Pixel> edge(pixels);
    Interpolate(a, b, edge);

    for (const Pixel &p : edge) {
      int y = p.y - minY;
      if (y < leftPixels.size()) {
        if (p.x < leftPixels[y].x) {
          leftPixels[y].x = p.x;
          leftPixels[y].zinv = p.zinv;
          leftPixels[y].pos3d = p.pos3d;
        }
        if (p.x > rightPixels[y].x) {
          rightPixels[y].x = p.x;
          rightPixels[y].zinv = p.zinv;
          rightPixels[y].pos3d = p.pos3d;
        }
      }
    }
  }
}

void DrawPolygonRows(const vector<Pixel> &leftPixels,
                     const vector<Pixel> &rightPixels) {
  for (size_t i = 0; i < leftPixels.size(); ++i) {
    if (leftPixels[i].x < rightPixels[i].x) {
      vector<Pixel> row(rightPixels[i].x - leftPixels[i].x + 1);
      Interpolate(leftPixels[i], rightPixels[i], row);
      for (const Pixel &p : row) {
        PixelShader(p);
      }
    }
  }
}

void DrawPolygon(const vector<Vertex> &vertices) {
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);
  for (int i = 0; i < V; ++i) {
    VertexShader(vertices[i], vertexPixels[i]);
  }
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(leftPixels, rightPixels);
}

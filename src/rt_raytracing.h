#pragma once

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>

namespace rt {

struct RTContext {
    int width = 500;
    int height = 500;
    std::vector<glm::vec4> image;
    bool freeze = false;
    int current_frame = 0;
    int current_line = 0;
    int max_frames = 1000;
    int max_bounces = 10;
    float epsilon = 2e-4f;
    glm::mat4 view = glm::mat4(1.0f);
    glm::vec3 ground_color = glm::vec3(1.0f, 1.0f, 1.0f);
    glm::vec3 sky_color = glm::vec3(0.5f, 0.7f, 1.0f);
    bool show_normals = false;
    // Add more settings and parameters here
    // ...
    bool gamma = true;
    bool anti_alias = true;
    glm::vec3 diffuseBallColor = glm::vec3(0.7f, 0.3f, 0.3f);
    glm::vec3 metalBallColor = glm::vec3(0.8f, 0.6f, 0.2f);
    glm::vec3 diffuseMeshColor = glm::vec3(0.5f, 0.75f, 0.5f);
};

void setupScene(RTContext &rtx, const char *mesh_filename);
void updateImage(RTContext &rtx);
void resetImage(RTContext &rtx);
void resetAccumulation(RTContext &rtx);

}  // namespace rt

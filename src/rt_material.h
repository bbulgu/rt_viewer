#pragma once

#include <glm/glm.hpp>
#include "rt_ray.h"
#include "rt_hitable.h"

namespace rt {
    struct hit_record;

    inline double random_double() {
        // Returns a random real in [0,1).
        return rand() / (RAND_MAX + 1.0);
    }

    inline double random_double(double min, double max) {
        // Returns a random real in [min,max).
        return min + (max - min) * random_double();
    }

    inline static glm::vec3 random() {
        return glm::vec3(random_double(), random_double(), random_double());
    }

    inline static glm::vec3 random_unit_vector() {
        return glm::normalize(random());
    }

    inline static glm::vec3 random(double min, double max) {
        return glm::vec3(random_double(min, max), random_double(min, max), random_double(min, max));
    }

    inline static double vector_length(glm::vec3 v)
    {
        return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    inline static glm::vec3 random_in_unit_sphere() {
        while (true) {
            glm::vec3 p = random(-1, 1);
            if (vector_length(p) * vector_length(p) >= 1) continue;
            return p;
        }
    }

    inline bool near_zero(glm::vec3 v) {
    // Return true if the vector is close to zero in all dimensions.
        const auto s = 1e-8;
        return (fabs(v[0]) < s) && (fabs(v[1]) < s) && (fabs(v[2]) < s);
    }
 

    class Material {
    public:
        virtual bool scatter(
            const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered
        ) const = 0;
    };

    class diffuse : public Material {
    public:
        diffuse(glm::vec3& a) : albedo(a) {}
        virtual bool scatter(
            const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered
        ) const override {
            glm::vec3 scatter_direction = rec.normal + random_unit_vector();

            // Catch degenerate scatter direction
            if (near_zero(scatter_direction))
                scatter_direction = rec.normal;

            scattered = Ray(rec.p, scatter_direction);
            attenuation = albedo;

            return true;
        }

    public:
        glm::vec3 albedo;
    };

    class metal : public Material {
    public:
        metal(const glm::vec3& a) : albedo(a) {}

        virtual bool scatter(
            const Ray& r_in, const HitRecord& rec, glm::vec3& attenuation, Ray& scattered
        ) const override {
            glm::vec3 reflected = glm::reflect(glm::normalize(r_in.direction()), rec.normal);
            scattered = Ray(rec.p, reflected);
            attenuation = albedo;
            return (glm::dot(scattered.direction(), rec.normal) > 0);
            /*
            glm::vec3 reflected = reflect(unit_vector(ray_in.direction()), rec.normal);
        scattered = Ray(rec.p, reflected);
        attenuation = albedo;

        return (glm::dot(scattered.direction(), rec.normal) > 0);
        */
        }

    public:
        glm::vec3 albedo;
    };

}
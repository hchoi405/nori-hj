/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
    /**
     * @brief A node of octree
     *
     */
    struct Node {
        Node *child[8];
        std::vector<uint32_t> triangles;
        BoundingBox3f bbox;

        Node(const BoundingBox3f &box) : bbox(box) {
            for (int i=0; i<8; ++i) child[i] = nullptr;
        }
    };

public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its,
                      bool shadowRay) const;

    /// Build octree using bounding box and triangles
    Node *buildOctree(const BoundingBox3f &bbox,
                      std::vector<uint32_t> &triangles);

    // Traverse the node with a ray
    bool traverseOctree(Node *node, Ray3f &ray) const;

private:
    Node *root;              ///< Root node of accel
    Mesh *m_mesh = nullptr;  ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;    ///< Bounding box of the entire scene
};

NORI_NAMESPACE_END

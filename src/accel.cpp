/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/accel.h>

#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh) throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() { /* Nothing to do here for now */
    std::vector<uint32_t> triangles;
    for (uint32_t t = 0; t < m_mesh->getTriangleCount(); ++t)
        triangles.push_back(t);
    root = buildOctree(m_mesh->getBoundingBox(), triangles);

    for (int i = 0; i < 8; ++i) {
    }
}

Accel::Node *Accel::buildOctree(const BoundingBox3f &bbox,
                                std::vector<uint32_t> &triangles) {
    if (triangles.empty()) return nullptr;

    if (triangles.size() < 10) {
        Node *node = new Node(bbox);
        node->triangles = triangles;
        return node;
    }

    std::vector<uint32_t> subTriangles[8];
    BoundingBox3f subBox[8];

    // Find bounding box for child
    for (int i = 0; i < 8; ++i) {
        subBox[i].min = bbox.getCenter().cwiseMin(bbox.getCorner(i));
        subBox[i].max = bbox.getCenter().cwiseMax(bbox.getCorner(i));
    }

    // Allocate triangles to child
    for (uint32_t t : triangles) {
        for (int i = 0; i < 8; ++i) {
            if (m_mesh->getBoundingBox(t).overlaps(subBox[i]))
                subTriangles[i].push_back(t);
        }
    }

    Node *node = new Node(bbox);
    for (int i = 0; i < 8; ++i) {
        node->child[i] = buildOctree(subBox[i], subTriangles[i]);
    }

    return node;
}

bool Accel::traverseOctree(Node *node, Ray3f &ray) const {
    if (node == nullptr) {
        return false;
    } else if (node->triangles.empty()) {
        return false;
    } else if (node->triangles.size() < 10) {
        return true;
    }

    bool intersect = false;
    if (node->bbox.rayIntersect(ray)) {
        for (int i = 0; i < 8; ++i) {
            if (node->child[i] != nullptr) {
                intersect = traverseOctree(node->child[i], ray);
            }
        }
    }

    return intersect;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its,
                         bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t)-1;  // Triangle index of the closest intersection

    Ray3f ray(ray_);  /// Make a copy of the ray (we will need to update its
    /// '.maxt' value)

    foundIntersection = traverseOctree(root, ray);

    // /* Brute force search through all triangles */
    // for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //     float u, v, t;
    //     if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //         /* An intersection was found! Can terminate
    //            immediately if this is a shadow ray query */
    //         if (shadowRay) return true;
    //         ray.maxt = its.t = t;
    //         its.uv = Point2f(u, v);
    //         its.mesh = m_mesh;
    //         f = idx;
    //         foundIntersection = true;
    //     }
    // }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1 - its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh = its.mesh;
        const MatrixXf &V = mesh->getVertexPositions();
        const MatrixXf &N = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) + bary.y() * UV.col(idx1) +
                     bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame =
                Frame((bary.x() * N.col(idx0) + bary.y() * N.col(idx1) +
                       bary.z() * N.col(idx2))
                          .normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

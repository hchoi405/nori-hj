/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/accel.h>
#include <nori/bbox.h>
#include <nori/mesh.h>
#include <stdio.h>

#include <Eigen/Geometry>
#include <iostream>
#include <list>
#include <vector>

using namespace std;

NORI_NAMESPACE_BEGIN

class Node {
 public:
  Node *childNodes[8] = {NULL};
  std::vector<uint32_t> triangles;
  BoundingBox3f bbox;

  Node() {}

  void addChildNode(Node *node, int i) {
    if (node != NULL) {
      childNodes[i] = node;
    }
  }

  void setTriangle(std::vector<uint32_t> triangleList) {
    // for(uint32_t i = 0; i < triangleList.size(); i++) {
    //   triangles.push_back(triangleList[i]);
    // }

    triangles = triangleList;
  }

  // std::vector<Node> getChildNodes () {
  //   return childNodes;
  // }

  bool isLeafNode() {
    if (childNodes[0] == NULL) {
      return true;
    } else {
      return false;
    }
  }
};

Node *buildOctree(BoundingBox3f &, std::vector<uint32_t> &, int, Mesh *);
// std::vector<uint32_t> triangleList;

class Tree {
 public:
  Node *root;
  std::vector<uint32_t> triangleList;

  Tree(Mesh *m_mesh) {
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); idx++) {
      triangleList.push_back(idx);
    }
    // root = buildOctree(m_bbox, triangleList, 1, m_mesh);
  }
};

Tree *tree;

void Accel::addMesh(Mesh *mesh) {
  if (m_mesh) throw NoriException("Accel: only a single mesh is supported!");
  m_mesh = mesh;
  m_bbox = m_mesh->getBoundingBox();
  tree = new Tree(m_mesh);
}

void Accel::build() {
  /* Nothing to do here for now */
  // std::vector<uint32_t> triangleList;

  // uint32_t idx = 0;
  // for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); idx++)
  // {
  //   triangleList.push_back(idx);
  // }
  // //printf("**%d\n", (int)triangleList.size());
  tree->root = buildOctree(m_bbox, tree->triangleList, 1, m_mesh);

  // tree = new Tree(m_mesh, m_bbox);
}

Node *buildOctree(BoundingBox3f &boundingBox,
                  std::vector<uint32_t> &triangleList, int depth,
                  Mesh *m_mesh) {
  Node *node = new Node();
  // printf("alloc done - %d\n", depth);

  // node -> setTriangle(triangleList);

  if (triangleList.empty()) {
    return node;  // null -> node
  } else if (triangleList.size() <= 10 || depth >= 3) {
    node->setTriangle(triangleList);
    node->bbox = boundingBox;
    // printf("------%d\n", (int)triangleList.size());
    return node;
  }

  std::vector<uint32_t> triList[8];
  BoundingBox3f subbox;
  // node->setTriangle(triangleList);
  node->bbox = boundingBox;

  for (int i = 0; i < 8; i++) {
    Vector3f minpoint(
        std::min(boundingBox.getCenter()[0], boundingBox.getCorner(i)[0]),
        std::min(boundingBox.getCenter()[1], boundingBox.getCorner(i)[1]),
        std::min(boundingBox.getCenter()[2], boundingBox.getCorner(i)[2]));
    Vector3f maxpoint(
        std::max(boundingBox.getCenter()[0], boundingBox.getCorner(i)[0]),
        std::max(boundingBox.getCenter()[1], boundingBox.getCorner(i)[1]),
        std::max(boundingBox.getCenter()[2], boundingBox.getCorner(i)[2]));
    subbox = BoundingBox3f(minpoint, maxpoint);

    for (uint32_t idx = 0; idx < triangleList.size(); idx++) {
      uint32_t triangle = triangleList[idx];
      // printf("1\n");
      BoundingBox3f triBound = m_mesh->getBoundingBox(triangle);
      // printf("2\n");

      if (triBound.overlaps(subbox, true)) {
        triList[i].push_back(triangle);
      }
    }

    node->addChildNode(buildOctree(subbox, triList[i], depth + 1, m_mesh), i);
  }

  return node;
}

bool traverse(const Ray3f ray_, Intersection &its, bool shadowRay, Node *curr,
              uint32_t &f, bool &foundIntersection, Mesh *mesh_, int depth) {
  Ray3f ray(ray_);
  // printf("traverse0\n");
  if (curr != NULL) {
    // printf("traverse1\n");
    if (!(*curr).isLeafNode()) {
      // printf("traverse2 %d\n", depth);
      for (int i = 0; i < 8; i++) {
        if ((*curr).childNodes[i] != nullptr)
          foundIntersection =
              traverse(ray_, its, shadowRay, (*curr).childNodes[i], f,
                       foundIntersection, mesh_, depth + 1);
      }
    } else {
      // printf("traverse3 %d\n", depth);
      // f = 0;
      //             ray.maxt = its.t = 1;
      // its.uv = Point2f(0, 0);
      // its.mesh = mesh_;
      // return true;
      float u, v, t;
      // printf("%d", (curr -> triangles).size());
      if ((curr->bbox).rayIntersect(ray)) {
        for (uint32_t i = 0; i < curr->triangles.size(); i++) {
          // printf("before triangle\n");
          uint32_t triangle = curr->triangles[i];
          // printf("after triangle\n");

          // printf("traverse4\n");
          if (mesh_->rayIntersect(
                  triangle, ray, u, v,
                  t)) {  // i -> triangle //compare with bbox? // bool
                         // rayIntersect(const Ray3f &ray) const
            // printf("traverse5\n");
            if (shadowRay) return true;

            if (ray.maxt > t) {
              ray.maxt = t;
              its.t = t;
              its.uv = Point2f(u, v);
              f = triangle;  // index
              its.mesh = mesh_;
              foundIntersection = true;
            }
          }
        }
      }
    }
  }

  return foundIntersection;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its,
                         bool shadowRay) const {
  bool foundIntersection = false;  // Was an intersection found so far?
  uint32_t f = (uint32_t)-1;       // Triangle index of the closest intersection

  Ray3f ray(ray_);  /// Make a copy of the ray (we will need to update its
                    /// '.maxt' value)

  /* Brute force search through all triangles */
  // for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
  //     float u, v, t;
  //     if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
  //         /* An intersection was found! Can terminate
  //            immediately if this is a shadow ray query */
  //         if (shadowRay)
  //             return true;
  //         ray.maxt = its.t = t;
  //         its.uv = Point2f(u, v);
  //         its.mesh = m_mesh;
  //         f = idx;
  //         foundIntersection = true;
  //     }
  // }

  foundIntersection = traverse(ray_, its, shadowRay, tree->root, f,
                               foundIntersection, m_mesh, 0);
  // printf("true false OK\n");

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

      its.shFrame = Frame((bary.x() * N.col(idx0) + bary.y() * N.col(idx1) +
                           bary.z() * N.col(idx2))
                              .normalized());
    } else {
      its.shFrame = its.geoFrame;
    }
  }

  return foundIntersection;
}

NORI_NAMESPACE_END

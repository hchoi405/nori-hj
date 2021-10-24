/*
        This file is part of Nori, a simple educational ray tracer

        Copyright (c) 2015 by Wenzel Jakob

        Nori is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License Version 3
        as published by the Free Software Foundation.

        Nori is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/frame.h>
#include <nori/vector.h>
#include <nori/warp.h>

#include <algorithm>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) { return sample; }

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
  return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f
                                                                      : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
  // throw NoriException("Warp::squareToTent() is not yet implemented!");
  float x = sample[0], y = sample[1];
  float px = (x < 0.5) ? (-2 + std::sqrt(8 * sample[0])) / 2
                       : px = (2 - std::sqrt(8 * (1 - sample[0]))) / 2;

  float py = (y < 0.5) ? (-2 + std::sqrt(8 * sample[1])) / 2
                       : (2 - std::sqrt(8 * (1 - sample[1]))) / 2;

  return Point2f(px, py);
}

float Warp::squareToTentPdf(const Point2f &p) {
  // throw NoriException("Warp::squareToTentPdf() is not yet implemented!");
  return (1 - std::abs(p[0])) * (1 - std::abs(p[1]));
}

float uniformToWeighted(float p) {
  // normalize
  p /= M_PI * 2;

  const float xrange[] = {0, 0.25f, 0.5f, 0.75f, 1.f};
  const float yrange[] = {0, 0.1f, 0.3f, 0.6f, 1.f};

  for (int i = 1; i < 5; ++i) {
    if (p < yrange[i]) {
      float ratio = (p - yrange[i - 1]) / (yrange[i] - yrange[i - 1]);
      return (ratio * (xrange[i] - xrange[i - 1]) + xrange[i - 1]) * M_PI * 2;
    }
  }
  return 0.f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
  // throw NoriException("Warp::squareToUniformDisk() is not yet implemented!");
  if (sample[0] == 0 && sample[1] == 0) {
    return Point2f(0, 0);
  }

  float r = 1;
  float theta = 0;

  float a = 2 * sample[0] - 1;
  float b = 2 * sample[1] - 1;

  if (a * a > b * b) {
    r *= a;
    theta = (M_PI / 4) * (b / a);
  } else {
    r *= b;
    theta = (M_PI / 2) - ((M_PI / 4) * (a / b));
  }

  float newX = r * std::cos(theta);
  float newY = r * std::sin(theta);

  theta = uniformToWeighted(std::atan2(newY, newX) + M_PI);

  newX = std::abs(r) * std::cos(theta);
  newY = std::abs(r) * std::sin(theta);

  return Point2f(newX, newY);
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
  // throw NoriException("Warp::squareToUniformDiskPdf() is not yet
  // implemented!");
  if (std::sqrt(p[0] * p[0] + p[1] * p[1]) <= 1)
    return 1 / M_PI;
  else
    return 0;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
  // throw NoriException("Warp::squareToUniformSphere() is not yet
  // implemented!");
  float pi = sample[0] * M_PI * 2;
  float theta = acos(2 * sample[1] - 1);
  float u = cos(theta);
  float x = sqrt(1 - u * u) * cos(pi);
  float y = sqrt(1 - u * u) * sin(pi);
  float z = u;
  return Vector3f(x, y, z);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
  // throw NoriException("Warp::squareToUniformSpherePdf() is not yet
  // implemented!");
  if (std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) <= 1)
    return 1 / (4 * M_PI);
  else
    return 0;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
  // throw NoriException("Warp::squareToUniformHemisphere() is not yet
  // implemented!");

  float pi = sample[0] * M_PI * 2;
  float theta = acos(sample[1]);
  float u = sample[1];
  float x = sqrt(1 - u * u) * cos(pi);
  float y = sqrt(1 - u * u) * sin(pi);
  float z = u;
  return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
  // throw NoriException("Warp::squareToUniformHemispherePdf() is not yet
  // implemented!");
  if (v[2] >= 0 && std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) <= 1)
    return 1 / (2 * M_PI);
  else
    return 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
  // throw NoriException("Warp::squareToCosineHemisphere() is not yet
  // implemented!");
  Point2f d = squareToUniformDisk(sample);
  float value = 1.f - d[0] * d[0] - d[1] * d[1];
  float z = std::sqrt(std::max(0.0f, value));
  return Vector3f(d[0], d[1], z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
  // throw NoriException("Warp::squareToCosineHemispherePdf() is not yet
  // implemented!");
  if (v[2] >= 0) {
    // printf("%f\n", acos(v.dot(Vector3f(0, 0, 1)) / v.norm()) * 180 / M_PI);
    return (v.dot(Vector3f(0, 0, 1)) / v.norm()) * INV_PI;
  } else
    return 0;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
  // throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
  float logsample = log(1 - sample.x());
  if (isinf(logsample)) {
    logsample = 0;
  }
  float tanSquare = -1 * alpha * alpha * logsample;
  float theta = atan(sqrt(tanSquare));
  float phi = 2 * M_PI * sample.y();
  return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

// float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
// 	// throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");

// 	float cosTheta = Frame::cosTheta(m);
// 	if (cosTheta <= 0) return 0;
// 	float sinTheta = Frame::sinTheta(m);
// 	float tanTheta = Frame::tanTheta(m);
// 	float tanSquare = tanTheta * tanTheta;
// 	//printf("%s\n", m.toString().c_str());

// 	// return (2 * exp(-pow(tanTheta, 2.f) / pow(alpha, 2.f))) / (pow(alpha, 2.f) * pow(cosTheta, 3.f)) / (2 * M_PI);
//     // return 2* sinTheta * exp((-1 * tanSquare)/ (alpha * alpha)) / (alpha * alpha * m.z() * m.z() * m.z());
// }

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
// throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
    float tanSquare = 1/(m.z() * m.z()) - 1.0f;
    float sinTheta = sqrt(1 - m.z() * m.z());
    if (m.z() > 0)
	    return (2 * exp(-tanSquare / pow(alpha, 2.f))) / (pow(alpha, 2.f) * pow(m.z(), 3.f)) / (2 * M_PI);
        // return 2* sinTheta * exp((-1 * tanSquare)/ (alpha * alpha)) / (alpha * alpha * m.z() * m.z() * m.z());
    else
        return 0;
}

NORI_NAMESPACE_END

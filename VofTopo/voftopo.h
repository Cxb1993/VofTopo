#ifndef VOFTOPO_H
#define VOFTOPO_H

#include <vector>
#include <cstdlib>
#include <utility>
#include <set>
#include <cmath>
#include "helper_math.h"

typedef int id_type;

typedef struct {
  float x;
  float y;
  float z;
  id_type id;
} f3u1_t;

inline 
f3u1_t make_f3u1_t(float x, float y, float z, unsigned int id)
{
  f3u1_t a = {x,y,z,id};
  return a;
}

//bool compF3u1_T(const f3u1_t &a, const f3u1_t &b);

template <typename T>
static void getGridPosition(f3u1_t particle, const int res[3],
			    const T *xcoords, 
			    const T *ycoords, 
			    const T *zcoords,
			    int idxOut[3], float bcoordOut[3])
{
  float3 origin = make_float3(xcoords[0], ycoords[0], zcoords[0]);
  float3 spacing = make_float3(xcoords[1] - xcoords[0], 
			       ycoords[1] - ycoords[0], 
			       zcoords[1] - zcoords[0]);
  float3 guessCoords = (make_float3(particle.x, particle.y, particle.z) - origin)/spacing;

  int idx[3] = {std::floor(guessCoords.x+0.5f),
		std::floor(guessCoords.y+0.5f),
		std::floor(guessCoords.z+0.5f)};
  const float prt[3] = {particle.x, particle.y, particle.z};
  const T *coords[3] = {xcoords, ycoords, zcoords};

  for (int c = 0; c < 3; c++) {

    while (idx[c] >= 0 && idx[c] < res[c]-1 && 
	   !(prt[c] >= coords[c][idx[c]] && prt[c] <= coords[c][idx[c]+1])) {

      if (prt[c] > coords[c][idx[c]+1]) {
	++idx[c];
      }
      else if (prt[c] < coords[c][idx[c]]) {
	--idx[c];
      }
    }

    if (idx[c] < 0) {
      idxOut[c] = 0;
      bcoordOut[c] = 0.0f;
    }
    else if (idx[c] >= res[c]-1) {
      idxOut[c] = res[c]-2;
      bcoordOut[c] = 1.0f;
    }
    else {
      idxOut[c] = idx[c];
      bcoordOut[c] = (prt[c] - coords[c][idx[c]])/(coords[c][idx[c]+1] - coords[c][idx[c]]);
    }
  }
}

template <typename T>
static float interpolateSca(const T *vofField,
			    const int* res, const int idxOut[3],
			    const float bcoordOut[3])
{
  int lx = idxOut[0];
  int ux = lx+1;
  int ly = idxOut[1];
  int uy = ly+1;
  int lz = idxOut[2];
  int uz = lz+1;

  float x = bcoordOut[0];
  float y = bcoordOut[1];
  float z = bcoordOut[2];

  unsigned lzslab = lz*res[0]*res[1];
  unsigned uzslab = uz*res[0]*res[1];
  int lyr = ly*res[0];
  int uyr = uy*res[0];

  unsigned id[8] = {lx + lyr + lzslab,
		    ux + lyr + lzslab,
		    lx + uyr + lzslab,
		    ux + uyr + lzslab,
		    lx + lyr + uzslab,
		    ux + lyr + uzslab,
		    lx + uyr + uzslab,
		    ux + uyr + uzslab};

  float vv[8];
  for (int i = 0; i < 8; i++) {
    vv[i] = vofField[id[i]];
  }

  float a = (1.0f-x)*vv[0] + x*vv[1];
  float b = (1.0f-x)*vv[2] + x*vv[3];
  float c = (1.0f-y)*a + y*b;
  a = (1.0f-x)*vv[4] + x*vv[5];
  b = (1.0f-x)*vv[6] + x*vv[7];
  float d = (1.0f-y)*a + y*b;

  return (1.0f-z)*c + z*d;
}

template <typename T>
static float3 interpolateVec(const T *velocityField,
			     const int* res, const int idxOut[3],
			     const float bcoordOut[3])
{
  int lx = idxOut[0];
  int ux = lx+1;
  int ly = idxOut[1];
  int uy = ly+1;
  int lz = idxOut[2];
  int uz = lz+1;

  float x = bcoordOut[0];
  float y = bcoordOut[1];
  float z = bcoordOut[2];

  unsigned lzslab = lz*res[0]*res[1];
  unsigned uzslab = uz*res[0]*res[1];
  int lyr = ly*res[0];
  int uyr = uy*res[0];

  unsigned id[8] = {lx + lyr + lzslab,
		    ux + lyr + lzslab,
		    lx + uyr + lzslab,
		    ux + uyr + lzslab,
		    lx + lyr + uzslab,
		    ux + lyr + uzslab,
		    lx + uyr + uzslab,
		    ux + uyr + uzslab};

  float3 vv[8];
  for (int i = 0; i < 8; i++) {
    vv[i].x = velocityField[id[i]*3+0];
    vv[i].y = velocityField[id[i]*3+1];
    vv[i].z = velocityField[id[i]*3+2];
  }

  float3 a = (1.0f-x)*vv[0] + x*vv[1];
  float3 b = (1.0f-x)*vv[2] + x*vv[3];
  float3 c = (1.0f-y)*a + y*b;
  a = (1.0f-x)*vv[4] + x*vv[5];
  b = (1.0f-x)*vv[6] + x*vv[7];
  float3 d = (1.0f-y)*a + y*b;

  return (1.0f-z)*c + z*d;
}

template <typename T>
void advectParticles(const T *velocityField, const int res[3], 
		     const T *xcoords, const T *ycoords, const T *zcoords, 
		     const float deltaT, std::vector<f3u1_t> &particles)
{
  std::vector<f3u1_t> particlesNext;
  std::vector<f3u1_t>::iterator it;
  for (it = particles.begin(); it != particles.end(); ++it) {

    f3u1_t particle = *it;
    if (particle.id > -1) {
      int idxOut[3];
      float bcoordOut[3];
      getGridPosition(particle, res, xcoords, ycoords, zcoords,
		      idxOut, bcoordOut);
      float3 velocity = interpolateVec(velocityField, res, idxOut, bcoordOut);
      float3 particleNext = make_float3(particle.x, particle.y, particle.z) + velocity*deltaT;
      particle = make_f3u1_t(particleNext.x, particleNext.y,
			     particleNext.z, it->id);
    }
    particlesNext.push_back(particle);
  }
  particles = particlesNext;
}

template <typename T>
void discardStrayParticles(const T *vofField, const int res[3], 
			   const T *xcoords, const T *ycoords, const T *zcoords, 
			   std::vector<f3u1_t> &particles)
{
  std::vector<f3u1_t> particlesValid;
  std::vector<f3u1_t>::iterator it;
  for (it = particles.begin(); it != particles.end(); ++it) {

    f3u1_t particle = *it;
    int idxOut[3];
    float bcoordOut[3];
    getGridPosition(particle, res, xcoords, ycoords, zcoords,
		    idxOut, bcoordOut);
    float f = interpolateSca(vofField, res, idxOut, bcoordOut);
    if (f <= 0.0f && particle.id > -1) {
      particle.id = particle.id*-1 - 1;
    }
    particlesValid.push_back(particle);
  }
  particles = particlesValid;  
}

#endif//VOFTOPO_H

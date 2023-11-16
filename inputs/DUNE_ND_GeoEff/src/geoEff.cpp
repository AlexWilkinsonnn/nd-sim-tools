#include "geoEff.h"
// C++ includes
#include <iostream>
#include <vector>
#include <random>
#include <math.h>

// Eigen Library
#include <Eigen/Dense>

geoEff::geoEff(int seed, bool verbose){

  verbosity = verbose;

  if (verbosity) {
    std::cout << "geoEff constructor called" << std::endl;
  }

  if (verbosity) {
    std::cout << "geoEff set seed to " << seed << std::endl;
  }

  N_THROWS = 64*64;
  N_THROWS_NDECC = 64*64;
  N_THROWS_FD = 64*64;

  if (verbosity){
    std::cout << "Number of throws at ND set to " << N_THROWS << std::endl;
    std::cout << "Number of throws at ND after ECC set to " << N_THROWS_NDECC << std::endl;
    std::cout << "Number of throws at FD set to " << N_THROWS_FD << std::endl;
  }

  prnGenerator = std::mt19937_64(seed);
  if (verbosity) {
    std::cout << "geoEff set random number generator to mt19937_64: seed= "<< seed << std::endl;
  }

  uniform = std::uniform_real_distribution<>(0., 1.);
  if (verbosity){
    std::cout << "geoEff set uniform distribution in [0, 1]" << std::endl;
  }

  translations[0].reserve(N_THROWS);
  translations[1].reserve(N_THROWS);
  translations[2].reserve(N_THROWS);
  rotations.reserve(N_THROWS);
  ndecctranslations[0].reserve(N_THROWS_NDECC);
  ndecctranslations[1].reserve(N_THROWS_NDECC);
  ndecctranslations[2].reserve(N_THROWS_NDECC);
  ndeccrotations.reserve(N_THROWS_NDECC);
  fdtranslations[0].reserve(N_THROWS_FD);
  fdtranslations[1].reserve(N_THROWS_FD);
  fdtranslations[2].reserve(N_THROWS_FD);

  if (verbosity){
    std::cout << "geoEff allocated memory for transformation vectors" << std::endl;
  }

  vetoSize = std::vector<float>(1,30.); // {30} cm
  vetoEnergy = std::vector<float>(1,30.); // 30 {MeV}

  if (verbosity){
    std::cout << "geoEff set veto size to 30 cm and energy threshold to 30 MeV" << std::endl;
  }

  useFixedBeamDir = false;

  // Initialize to all dimensions randomized
  for (int dim = 0; dim < 3; dim++) {
    randomizeVertex[dim] = true;
    randomizeVertexfd[dim] = true;
  }
}

void geoEff::setNthrows(unsigned long n){
  N_THROWS = n;
  if (verbosity){
    std::cout << "geoEff set number of throws at ND to " << N_THROWS << std::endl;
  }

  if (N_THROWS%64 && verbosity) std::cout << "geoEff warning: number of throws should be multiple of 64 for optimal use of output format."  << std::endl;
}

void geoEff::setNthrowsNDECC(unsigned long n){
  N_THROWS_NDECC = n;
  if (verbosity){
    std::cout << "geoEff set number of throws at FD to " << N_THROWS_NDECC << std::endl;
  }

  if (N_THROWS_NDECC%64 && verbosity) std::cout << "geoEff warning: number of throws should be multiple of 64 for optimal use of output format."  << std::endl;
}

void geoEff::setNthrowsFD(unsigned long n){
  N_THROWS_FD = n;
  if (verbosity){
    std::cout << "geoEff set number of throws at FD to " << N_THROWS_FD << std::endl;
  }

  if (N_THROWS_FD%64 && verbosity) std::cout << "geoEff warning: number of throws should be multiple of 64 for optimal use of output format."  << std::endl;
}

void geoEff::setVertex(float x, float y, float z){
  vertex[0] = x;
  vertex[1] = y;
  vertex[2] = z;
  if(verbosity){
    std::cout << "geoEff set vertex at ND to " << vertex[0] << " "<< vertex[1] << " "<< vertex[2] << std::endl;
  }
}

void geoEff::setNDrandVertex(float x, float y, float z){
  ndrandvertex[0] = x;
  ndrandvertex[1] = y;
  ndrandvertex[2] = z;
  if(verbosity){
    std::cout << "geoEff set nd event random throw vertex to " << ndrandvertex[0] << " "<< ndrandvertex[1] << " "<< ndrandvertex[2] << std::endl;
  }
}

void geoEff::setVertexFD(float x, float y, float z){
  fdvertex[0] = x;
  fdvertex[1] = y;
  fdvertex[2] = z;
  if(verbosity){
    std::cout << "geoEff set vertex at FD to " << fdvertex[0] << " "<< fdvertex[1] << " "<< fdvertex[2] << std::endl;
  }
}

void geoEff::setHitSegEdeps(std::vector<float> thishitSegEdeps){
  hitSegEdeps = thishitSegEdeps;
  if (verbosity) {
    std::cout << "geoEff setting hit segment energy deposits to ";
    for (unsigned int i = 0; i < hitSegEdeps.size(); i++) std::cout << hitSegEdeps[i] << " ";
    std::cout << std::endl;
  }
}

void geoEff::setHitSegPoss(std::vector<float> thishitSegPoss){

  // Set the vector
  hitSegPoss = thishitSegPoss;
  if (verbosity) {
    std::cout << "geoEff setting hit segment positions to ";
    for (unsigned int i = 0; i < hitSegPoss.size(); i++) std::cout << hitSegPoss[i] << " ";
    std::cout << std::endl;
  }

}

void geoEff::setRangeX(float xmin, float xmax){
  range[0][0] = xmin;
  range[0][1] = xmax;
}
void geoEff::setRangeY(float ymin, float ymax){
  range[1][0] = ymin;
  range[1][1] = ymax;
}
void geoEff::setRangeZ(float zmin, float zmax){
  range[2][0] = zmin;
  range[2][1] = zmax;
}

void geoEff::setRangeXFD(float xmin, float xmax){
  fdrange[0][0] = xmin;
  fdrange[0][1] = xmax;
}
void geoEff::setRangeYFD(float ymin, float ymax){
  fdrange[1][0] = ymin;
  fdrange[1][1] = ymax;
}
void geoEff::setRangeZFD(float zmin, float zmax){
  fdrange[2][0] = zmin;
  fdrange[2][1] = zmax;
}

void geoEff::setRandomizeX(bool r){
  randomizeVertex[0] = r;
}

void geoEff::setRandomizeY(bool r){
  randomizeVertex[1] = r;
}

void geoEff::setRandomizeZ(bool r){
  randomizeVertex[2] = r;
}

void geoEff::setRandomizeXFD(bool r){
  randomizeVertexfd[0] = r;
}

void geoEff::setRandomizeYFD(bool r){
  randomizeVertexfd[1] = r;
}

void geoEff::setRandomizeZFD(bool r){
  randomizeVertexfd[2] = r;
}

void geoEff::setActiveX(float xmin, float xmax){
  active[0][0] = xmin;
  active[0][1] = xmax;
}
void geoEff::setActiveY(float ymin, float ymax){
  active[1][0] = ymin;
  active[1][1] = ymax;
}
void geoEff::setActiveZ(float zmin, float zmax){
  active[2][0] = zmin;
  active[2][1] = zmax;
}

void geoEff::setFDActiveX(float xmin, float xmax){
  fdactive[0][0] = xmin;
  fdactive[0][1] = xmax;
}
void geoEff::setFDActiveY(float ymin, float ymax){
  fdactive[1][0] = ymin;
  fdactive[1][1] = ymax;
}
void geoEff::setFDActiveZ(float zmin, float zmax){
  fdactive[2][0] = zmin;
  fdactive[2][1] = zmax;
}

void geoEff::setOffsetX(float x){
  offset[0] = x;
}

void geoEff::setOffsetY(float y){
  offset[1] = y;
}

void geoEff::setOffsetZ(float z){
  offset[2] = z;
}

void geoEff::setBeamDir(float xdir, float ydir, float zdir){
  beamdir[0] = xdir;
  beamdir[1] = ydir;
  beamdir[2] = zdir;
}

void geoEff::setDecayPos(float x, float y, float z){
  decaypos[0] = x;
  decaypos[1] = y;
  decaypos[2] = z;
}

// For creating paired dataset for near to far det translation training
// Need a new decay position for random throws where vertex x is changed
void geoEff::setDecayPos4RandomThrowX(float x, float y, float z){
  decayposrandthrow[0] = x;
  decayposrandthrow[1] = y;
  decayposrandthrow[2] = z;
}

float geoEff::getDecayPos(int dim)
{
  return decaypos[dim];
}

void geoEff::setUseFixedBeamDir(bool use){
  useFixedBeamDir = use;
}

void geoEff::setVetoSizes(std::vector< float > vSizes){
  vetoSize = vSizes;
  if(verbosity){
    std::cout << "geoEff set veto sizes to ";
    for (unsigned int i = 0; i < vetoSize.size(); i++) std::cout << vetoSize[i] << " ";
  }
  std::cout << std::endl;
}
void geoEff::setVetoEnergyThresholds(std::vector< float > vThresholds){
  vetoEnergy = vThresholds;
  if(verbosity){
    std::cout << "geoEff set veto energy thresholds to ";
    for (unsigned int i = 0; i < vetoEnergy.size(); i++) std::cout << vetoEnergy[i] << " ";
  }
  std::cout << std::endl;
}

void geoEff::throwTransforms(){

  // Clear vectors
  translations[0].clear();
  translations[1].clear();
  translations[2].clear();
  rotations.clear();


  // dim from 0 to 2, corresponding to x, y and z
  for (int dim = 0; dim < 3; dim++){
    if (not randomizeVertex[dim]){
      translations[dim].resize(0,0);
    } else {
      translations[dim].clear();
      for (unsigned int i = 0; i < N_THROWS; i++){
        translations[dim].emplace_back(uniform(prnGenerator)*(range[dim][1]-range[dim][0])+range[dim][0]+offset[dim]);
      }
    }
  }

  rotations.clear();
  for (unsigned int i = 0; i < N_THROWS; i++){
    rotations.emplace_back((uniform(prnGenerator)-0.5)*2*M_PI);
  }

}

void geoEff::throwTransformsNDECC(){

  // Clear vectors
  ndecctranslations[0].clear();
  ndecctranslations[1].clear();
  ndecctranslations[2].clear();
  ndeccrotations.clear();

  // Rotate around on axis ND LAr center or FD beam axis, keep same orientation after ECC
  for (unsigned int i = 0; i < N_THROWS_NDECC; i++){
    ndeccrotations.emplace_back((uniform(prnGenerator)-0.5)*2*M_PI);
  }

  // dim from 0 to 2, corresponding to x, y and z
  for (int dim = 0; dim < 3; dim++){
    if (not randomizeVertexfd[dim]){
      ndecctranslations[dim].resize(0,0);
    } else {
      ndecctranslations[dim].clear();
      for (unsigned int i = 0; i < N_THROWS_NDECC; i++){
        ndecctranslations[dim].emplace_back(uniform(prnGenerator)*(range[dim][1]-range[dim][0])+range[dim][0]+offset[dim]);
      }
    }
  }

}

void geoEff::throwTransformsFD(){

  // Clear vectors
  fdtranslations[0].clear();
  fdtranslations[1].clear();
  fdtranslations[2].clear();

  // dim from 0 to 2, corresponding to x, y and z
  for (int dim = 0; dim < 3; dim++){
    if (not randomizeVertexfd[dim]){
      fdtranslations[dim].resize(0,0);
    } else {
      fdtranslations[dim].clear();
      for (unsigned int i = 0; i < N_THROWS_FD; i++){
        fdtranslations[dim].emplace_back(uniform(prnGenerator)*(fdrange[dim][1]-fdrange[dim][0])+fdrange[dim][0]);
      }
    }
  }

  // no need to rotate, keep same orientation from ND throw
}

std::vector< Eigen::Transform<float,3,Eigen::Affine> > geoEff::getTransforms(unsigned int iStart, int iEnd){

  unsigned int thisEnd;
  if (iEnd < 0) thisEnd = N_THROWS;
  else if (iEnd >= 0){
    thisEnd = iEnd;
    if (thisEnd > N_THROWS) thisEnd = N_THROWS;
  }

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms;

  // Tranformations that do not depend on the throws:
  // Move vertex to coordinate system origin to apply rotation
  Eigen::Affine3f tThere(Eigen::Translation3f(Eigen::Vector3f(-vertex[0], -vertex[1], -vertex[2])));
  // Move vertex back
  Eigen::Affine3f tBack(Eigen::Translation3f(Eigen::Vector3f(vertex[0], vertex[1], vertex[2])));

  for (unsigned int iThrow = iStart; iThrow < thisEnd; iThrow++){

    // Vertex displacement:
    Eigen::Affine3f tThrow(Eigen::Translation3f(Eigen::Vector3f(randomizeVertex[0] ? translations[0][iThrow]-vertex[0] : 0.,
								randomizeVertex[1] ? translations[1][iThrow]-vertex[1] : 0.,
								randomizeVertex[2] ? translations[2][iThrow]-vertex[2] : 0.)));

    // Rotation
    Eigen::Affine3f rThrow;
    if (useFixedBeamDir){
      rThrow = Eigen::Affine3f(Eigen::AngleAxisf(rotations[iThrow], Eigen::Vector3f(beamdir[0], beamdir[1], beamdir[2])));
    } else {
      // Calculate rotation due to translation
      // Calculate rotation angle
      double decayToVertex[3] = {0};
      double decayToTranslated[3] = {0};
      double translationAngle = 0, magDecayToVertex = 0, magDecayToTranslated = 0;
      for (int dim = 0; dim < 3; dim++) {
        decayToVertex[dim] = vertex[dim]-decaypos[dim];
        decayToTranslated[dim] = randomizeVertex[dim] ? translations[dim][iThrow]-decaypos[dim] : vertex[dim]-decaypos[dim];

        translationAngle += (decayToVertex[dim])*(decayToTranslated[dim]);
        magDecayToVertex += pow(decayToVertex[dim], 2);
        magDecayToTranslated += pow(decayToTranslated[dim], 2);
      }
      magDecayToVertex = sqrt(magDecayToVertex);
      magDecayToTranslated = sqrt(magDecayToTranslated);
      translationAngle /= (magDecayToVertex*magDecayToTranslated);
      translationAngle = acos(translationAngle);

      // Calculate rotation axis
      // Cross-product
      float translationAxis[3] = {0};
      translationAxis[0] = decayToVertex[1]*decayToTranslated[2] - decayToVertex[2]*decayToTranslated[1];
      translationAxis[1] = decayToVertex[2]*decayToTranslated[0] - decayToVertex[0]*decayToTranslated[2];
      translationAxis[2] = decayToVertex[0]*decayToTranslated[1] - decayToVertex[1]*decayToTranslated[0];
      float magTranslationAxis = 0.;
      for (int dim = 0; dim < 3; dim++) magTranslationAxis += pow(translationAxis[dim], 2);
      magTranslationAxis = sqrt(magTranslationAxis);

      if(magTranslationAxis!=0)  {for (int dim = 0; dim < 3; dim++) translationAxis[dim] /= magTranslationAxis;}
      else{for (int dim = 0; dim < 3; dim++) {translationAxis[dim] = 0.;translationAngle =0.;}}
      Eigen::Affine3f rTranslation(Eigen::Affine3f(Eigen::AngleAxisf(translationAngle, Eigen::Vector3f(translationAxis[0], translationAxis[1], translationAxis[2]))));

      // Calculate rotation due to thrown angle
      Eigen::Affine3f rPhiThrow(Eigen::Affine3f(Eigen::AngleAxisf(rotations[iThrow], Eigen::Vector3f(decayToTranslated[0]/magDecayToTranslated, decayToTranslated[1]/magDecayToTranslated, decayToTranslated[2]/magDecayToTranslated))));

      // Combine
      rThrow = rPhiThrow * rTranslation;
    }

    // Put everything together in single transform and store.
    transforms.emplace_back(tThrow * tBack * rThrow * tThere);
  }

  return transforms;
}

std::vector< Eigen::Transform<float,3,Eigen::Affine> > geoEff::getTransforms4RandomThrowX(unsigned int iStart, int iEnd){

  unsigned int thisEnd;
  if (iEnd < 0) thisEnd = N_THROWS;
  else if (iEnd >= 0){
    thisEnd = iEnd;
    if (thisEnd > N_THROWS) thisEnd = N_THROWS;
  }

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms;

  // Tranformations that do not depend on the throws:
  // Move vertex to coordinate system origin to apply rotation
  Eigen::Affine3f tThere(Eigen::Translation3f(Eigen::Vector3f(-vertex[0], -vertex[1], -vertex[2])));
  // Move vertex back
  Eigen::Affine3f tBack(Eigen::Translation3f(Eigen::Vector3f(vertex[0], vertex[1], vertex[2])));

  for (unsigned int iThrow = iStart; iThrow < thisEnd; iThrow++){

    // Vertex displacement:
    Eigen::Affine3f tThrow(Eigen::Translation3f(Eigen::Vector3f(randomizeVertex[0] ? translations[0][iThrow]-vertex[0] : 0.,
								randomizeVertex[1] ? translations[1][iThrow]-vertex[1] : 0.,
								randomizeVertex[2] ? translations[2][iThrow]-vertex[2] : 0.)));

    // Rotation
    Eigen::Affine3f rThrow;
    if (useFixedBeamDir){
      rThrow = Eigen::Affine3f(Eigen::AngleAxisf(rotations[iThrow], Eigen::Vector3f(beamdir[0], beamdir[1], beamdir[2])));
    } else {
      // Calculate rotation due to translation
      // Calculate rotation angle
      double decayToVertex[3] = {0};
      double decayToTranslated[3] = {0};
      double translationAngle = 0, magDecayToVertex = 0, magDecayToTranslated = 0;
      for (int dim = 0; dim < 3; dim++) {
        decayToVertex[dim] = vertex[dim]-decaypos[dim];
        decayToTranslated[dim] = randomizeVertex[dim] ? translations[dim][iThrow]-decayposrandthrow[dim] : vertex[dim]-decayposrandthrow[dim];

        translationAngle += (decayToVertex[dim])*(decayToTranslated[dim]);
        magDecayToVertex += pow(decayToVertex[dim], 2);
        magDecayToTranslated += pow(decayToTranslated[dim], 2);
      }
      magDecayToVertex = sqrt(magDecayToVertex);
      magDecayToTranslated = sqrt(magDecayToTranslated);
      translationAngle /= (magDecayToVertex*magDecayToTranslated);
      translationAngle = acos(translationAngle);

      // Calculate rotation axis
      // Cross-product
      float translationAxis[3] = {0};
      translationAxis[0] = decayToVertex[1]*decayToTranslated[2] - decayToVertex[2]*decayToTranslated[1];
      translationAxis[1] = decayToVertex[2]*decayToTranslated[0] - decayToVertex[0]*decayToTranslated[2];
      translationAxis[2] = decayToVertex[0]*decayToTranslated[1] - decayToVertex[1]*decayToTranslated[0];
      float magTranslationAxis = 0.;
      for (int dim = 0; dim < 3; dim++) magTranslationAxis += pow(translationAxis[dim], 2);
      magTranslationAxis = sqrt(magTranslationAxis);

      if(magTranslationAxis!=0)  {for (int dim = 0; dim < 3; dim++) translationAxis[dim] /= magTranslationAxis;}
      else{for (int dim = 0; dim < 3; dim++) {translationAxis[dim] = 0.;translationAngle =0.;}}
      Eigen::Affine3f rTranslation(Eigen::Affine3f(Eigen::AngleAxisf(translationAngle, Eigen::Vector3f(translationAxis[0], translationAxis[1], translationAxis[2]))));

      // Calculate rotation due to thrown angle
      Eigen::Affine3f rPhiThrow(Eigen::Affine3f(Eigen::AngleAxisf(rotations[iThrow], Eigen::Vector3f(decayToTranslated[0]/magDecayToTranslated, decayToTranslated[1]/magDecayToTranslated, decayToTranslated[2]/magDecayToTranslated))));

      // Combine
      rThrow = rPhiThrow * rTranslation;
    }

    // Put everything together in single transform and store.
    transforms.emplace_back(tThrow * tBack * rThrow * tThere);
  }

  return transforms;
}

std::vector< Eigen::Transform<float,3,Eigen::Affine> > geoEff::getTransformsFD(unsigned int iStart, int iEnd){

  unsigned int thisEnd;
  if (iEnd < 0) thisEnd = N_THROWS_FD;
  else if (iEnd >= 0){
    thisEnd = iEnd;
    if (thisEnd > N_THROWS_FD) thisEnd = N_THROWS_FD;
  }

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transformsfd;

  for (unsigned int iThrow = iStart; iThrow < thisEnd; iThrow++){

    // Vertex displacement:
    Eigen::Affine3f tThrow(Eigen::Translation3f(Eigen::Vector3f(randomizeVertexfd[0] ? fdtranslations[0][iThrow]-fdvertex[0] : 0.,
								randomizeVertexfd[1] ? fdtranslations[1][iThrow]-fdvertex[1] : 0.,
								randomizeVertexfd[2] ? fdtranslations[2][iThrow]-fdvertex[2] : 0.)));

    // No need to rotate, keep same orientation from ND rand throw

    // Put everything together in single transform and store.
    transformsfd.emplace_back(tThrow);
  }

  return transformsfd;
}

std::vector<float> geoEff::getCurrentThrowTranslationsX(){
  return translations[0];
}
std::vector<float> geoEff::getCurrentThrowTranslationsY(){
  return translations[1];
}
std::vector<float> geoEff::getCurrentThrowTranslationsZ(){
  return translations[2];
}
std::vector<float> geoEff::getCurrentThrowRotations(){
  return rotations;
}

std::vector<float> geoEff::getCurrentNDECCThrowTranslationsX(){
  return ndecctranslations[0];
}
std::vector<float> geoEff::getCurrentNDECCThrowTranslationsY(){
  return ndecctranslations[1];
}
std::vector<float> geoEff::getCurrentNDECCThrowTranslationsZ(){
  return ndecctranslations[2];
}
std::vector<float> geoEff::getCurrentNDECCThrowRotations(){
  return ndeccrotations;
}

std::vector<float> geoEff::getCurrentFDThrowTranslationsX(){
  return fdtranslations[0];
}
std::vector<float> geoEff::getCurrentFDThrowTranslationsY(){
  return fdtranslations[1];
}
std::vector<float> geoEff::getCurrentFDThrowTranslationsZ(){
  return fdtranslations[2];
}


// Get the coordinates of hadron hits after eigen transformation, i is the # of throw
std::vector< float > geoEff::getCurrentThrowDeps(int i, int dim){

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));

  Eigen::Matrix3Xf transformedEdeps = getTransforms(i, i+1)[0] * hitSegPosOrig;

  int nEdeps = hitSegEdeps.size();

  std::vector< float > ret(nEdeps);

  for (int iDep = 0; iDep < nEdeps; iDep++){
    ret[iDep] = transformedEdeps(dim, iDep);
  }

  return ret;
}

std::vector< float > geoEff::getCurrentThrowDepsX(int i){
  return getCurrentThrowDeps(i, 0);
}
std::vector< float > geoEff::getCurrentThrowDepsY(int i){
  return getCurrentThrowDeps(i, 1);
}
std::vector< float > geoEff::getCurrentThrowDepsZ(int i){
  return getCurrentThrowDeps(i, 2);
}


std::vector< std::vector< std::vector< uint64_t > > > geoEff::getHadronContainmentThrows(bool ignore_uncontained){

  // Figure out how many multiples of 64 bits needed to store output
  int n_longs = N_THROWS / 64;
  if (N_THROWS % 64) n_longs++;

  // Pass/fail for each set of vetoSize and vetoEnergy
  std::vector< std::vector< std::vector< uint64_t > > > hadronContainment(vetoSize.size(), std::vector< std::vector< uint64_t > >(vetoEnergy.size(), std::vector < uint64_t >(n_longs, 0)));

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));

  // Check if event is contained by any of the existing conditions
  if (ignore_uncontained) {
    int origContained = 0;
    std::vector< std::vector< bool > > vecOrigContained = getHadronContainmentOrigin();
    for (unsigned int i = 0; i < vetoSize.size(); i++){
      for (unsigned int j = 0; j < vetoEnergy.size(); j++){
        if (vecOrigContained[i][j]) origContained++;
      }
    }

    // If not, then return
    if (origContained == 0) return hadronContainment;
  }

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms = getTransforms();
  // Else, loop through set of rotation translations
  for (unsigned int t = 0; t < N_THROWS; t++){
    // Apply transformation to energy deposit positions
    Eigen::Matrix3Xf transformedEdeps = transforms[t] * hitSegPosOrig;
    // Loop through conditions
    for (unsigned int i = 0; i < vetoSize.size(); i++){
      for (unsigned int j = 0; j < vetoEnergy.size(); j++){
        // Check containment and set bit
        if (isContained(transformedEdeps, hitSegEdeps, vetoSize[i], vetoEnergy[j]))
        {
	         hadronContainment[i][j][t/64] |= ((uint64_t)1)<<(t%64);
	      }
      }
    }
  }

  return hadronContainment;
}

struct throwcombo geoEff::getNDContainment4RandomThrowX(){

  struct throwcombo ndthrowcombo;

  // Figure out how many multiples of 64 bits needed to store output
  int n_longs = N_THROWS / 64;
  if (N_THROWS % 64) n_longs++;

  // Pass/fail for each set of vetoSize and vetoEnergy
  std::vector< std::vector< std::vector< uint64_t > > > NDContainment4RandomThrowX(vetoSize.size(), std::vector< std::vector< uint64_t > >(vetoEnergy.size(), std::vector < uint64_t >(n_longs, 0)));

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));
  std::vector< Eigen::Matrix3Xf > transformedEdepss;

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms4RandomThrowX = getTransforms4RandomThrowX();
  // Else, loop through set of rotation translations
  for (unsigned int t = 0; t < N_THROWS; t++){
    // Apply transformation to energy deposit positions
    Eigen::Matrix3Xf transformedEdeps4RandomThrowX = transforms4RandomThrowX[t] * hitSegPosOrig;
    transformedEdepss.emplace_back(transformedEdeps4RandomThrowX);
    // Loop through conditions
    for (unsigned int i = 0; i < vetoSize.size(); i++){
      for (unsigned int j = 0; j < vetoEnergy.size(); j++){
        // Check containment and set bit
        if (isContainedInND(transformedEdeps4RandomThrowX, hitSegEdeps, vetoSize[i], vetoEnergy[j]))
        {
	         NDContainment4RandomThrowX[i][j][t/64] |= ((uint64_t)1)<<(t%64);
	      }
      }
    }
  }

  ndthrowcombo.thrownEdepspos = transformedEdepss;
  ndthrowcombo.containresult = NDContainment4RandomThrowX;

  return ndthrowcombo;
}

struct throwcombo geoEff::getFDContainment4RandomThrow(Eigen::Matrix3Xf FDhitSegPosOrig){

  struct throwcombo fdthrowcombo;

  // Figure out how many multiples of 64 bits needed to store output
  int n_longs = N_THROWS_FD / 64;
  if (N_THROWS_FD % 64) n_longs++;

  std::vector< Eigen::Matrix3Xf > transformedEdepssFD;

  // Pass/fail for each set of vetoSize and vetoEnergy
  std::vector< std::vector< std::vector< uint64_t > > > FDContainment4RandomThrow(vetoSize.size(), std::vector< std::vector< uint64_t > >(vetoEnergy.size(), std::vector < uint64_t >(n_longs, 0)));

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transformsFD4RandomThrow = getTransformsFD();
  // Else, loop through set of rotation translations
  for (unsigned int t = 0; t < N_THROWS_FD; t++){
    // Apply transformation to energy deposit positions
    Eigen::Matrix3Xf transformedFDEdeps4RandomThrow = transformsFD4RandomThrow[t] * FDhitSegPosOrig;
    transformedEdepssFD.emplace_back(transformedFDEdeps4RandomThrow);
    // Loop through conditions
    for (unsigned int i = 0; i < vetoSize.size(); i++){
      for (unsigned int j = 0; j < vetoEnergy.size(); j++){
        // Check containment and set bit
        if (isPairedFDEvtContainedInFD(transformedFDEdeps4RandomThrow, hitSegEdeps, vetoSize[i], vetoEnergy[j]))
        {
	         FDContainment4RandomThrow[i][j][t/64] |= ((uint64_t)1)<<(t%64);
	      }
      }
    }
  }

  fdthrowcombo.thrownEdepspos = transformedEdepssFD;
  fdthrowcombo.containresult = FDContainment4RandomThrow;

  return fdthrowcombo;
}

void geoEff::setSeed(int seed){
  prnGenerator = std::mt19937_64(seed);
}

std::vector< std::vector< bool > > geoEff::getHadronContainmentOrigin(){
  // Initialize return vector
  std::vector< std::vector< bool > > hadronContainment(vetoSize.size(), std::vector< bool >(vetoEnergy.size(), false));

  // Set Eigen Map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > hitSegPosOrig(hitSegPoss.data(),3,hitSegPoss.size()/3,Eigen::OuterStride<>(3));

  for (unsigned int i = 0; i < vetoSize.size(); i++){
    for (unsigned int j = 0; j < vetoEnergy.size(); j++){
      if (isContained(hitSegPosOrig, hitSegEdeps, vetoSize[i], vetoEnergy[j])) hadronContainment[i][j] = true;
    }
  }

  return hadronContainment;
}

void geoEff::setOffAxisOffsetX(float x){
  OffAxisOffset[0] = x;
}

void geoEff::setOffAxisOffsetY(float y){
  OffAxisOffset[1] = y;
}

void geoEff::setOffAxisOffsetZ(float z){
  OffAxisOffset[2] = z;
}

bool geoEff::isContained( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits, float vSize, float vetoEnergyThreshold ){

  float vetoEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      // low
      if ( (hitSegments(dim, i)-OffAxisOffset[dim] < active[dim][0]+vSize) and
           (hitSegments(dim, i)-OffAxisOffset[dim] > active[dim][0]) ) {
        vetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
      // high
      if ( (hitSegments(dim, i)-OffAxisOffset[dim] > active[dim][1]-vSize) and
           (hitSegments(dim, i)-OffAxisOffset[dim] < active[dim][1]) ) {
        vetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
    }
  }

  return vetoEnergy < vetoEnergyThreshold;
}

// New function added for N2FD study
bool geoEff::isContainedInND( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits, float vSize, float vetoEnergyThreshold ){

  float vetoEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      // low
      if ( (hitSegments(dim, i)-offset[dim] < active[dim][0]+vSize) and
           (hitSegments(dim, i)-offset[dim] > active[dim][0]) ) {
        vetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
      // high
      if ( (hitSegments(dim, i)-offset[dim] > active[dim][1]-vSize) and
           (hitSegments(dim, i)-offset[dim] < active[dim][1]) ) {
        vetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
    }
  }

  return vetoEnergy < vetoEnergyThreshold;
}

// New function added for N2FD study
bool geoEff::isPairedFDEvtContainedInFD( Eigen::Matrix3Xf hitSegments, std::vector<float> energyDeposits, float vSize, float vetoEnergyThreshold ){

  float fdvetoEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    for (int dim = 0; dim < 3; dim++){
      // low
      if ( (hitSegments(dim, i) < fdactive[dim][0]+vSize) and
           (hitSegments(dim, i) > fdactive[dim][0]) ) {
        fdvetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
      // high
      if ( (hitSegments(dim, i) > fdactive[dim][1]-vSize) and
           (hitSegments(dim, i) < fdactive[dim][1]) ) {
        fdvetoEnergy += energyDeposits[i];
        break; // Only count each energy deposit once
      }
    }
  }

  return fdvetoEnergy < vetoEnergyThreshold;
}

// Get TOTAL E
float geoEff::getTotE(std::vector<float> energyDeposits){

  float totEnergy = 0.;

  for (unsigned int i = 0; i < energyDeposits.size(); i++){
    totEnergy += energyDeposits[i];
  }
  return totEnergy;
}

float geoEff::getCurrentThrowsTotE(){
  float totEnergy = 0.;
  totEnergy = getTotE(hitSegEdeps);
  return totEnergy;
}


//
// #########################################################################
// #########################################################################
// #########################################################################
// #########################################################################
// #########################################################################
// #########################################################################
// #########################################################################
// #########################################################################
// Translation from On-Axis position to Off-Axis positions
void geoEff::setOnAxisVertex(float x, float y, float z){
  OnAxisVertex[0] = x;
  OnAxisVertex[1] = y;
  OnAxisVertex[2] = z;
  if(verbosity){
    std::cout << "geoEff set On-Axis vertex to " << OnAxisVertex[0] << " "<< OnAxisVertex[1] << " "<< OnAxisVertex[2] << std::endl;
  }
}
// Vertex before rotations
void geoEff::setOffAxisVertex(float x, float y, float z){
  OffAxisVertex[0] = x;
  OffAxisVertex[1] = y;
  OffAxisVertex[2] = z;
}


// Position vector space Eigen transformation
std::vector< Eigen::Transform<float,3,Eigen::Affine> > geoEff::getTransforms_NDtoND(){

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms_NDtoND;

  // Tranformations that do not depend on the throws:
  // Move vertex to coordinate system origin to apply rotation
  Eigen::Affine3f tThere_NDtoND(Eigen::Translation3f(Eigen::Vector3f(-OffAxisVertex[0], -OffAxisVertex[1], -OffAxisVertex[2])));
  // Move vertex back
  Eigen::Affine3f tBack_NDtoND(Eigen::Translation3f(Eigen::Vector3f(OffAxisVertex[0], OffAxisVertex[1], OffAxisVertex[2])));
  // Eigen::Affine3f is a typedef of Eigen::Transform<float, 3, Eigen::Affine>

    // Vertex displacement:
    // Eigen::Affine3f tThrow_NDtoND(Eigen::Translation3f(Eigen::Vector3f(OffAxisVertex[0]-OnAxisVertex[0],OffAxisVertex[1]-OnAxisVertex[1],OffAxisVertex[2]-OnAxisVertex[2])));

    // Rotation
    Eigen::Affine3f rThrow_NDtoND;
    {
      // Calculate rotation due to translation
      // Calculate rotation angle
      double decayToVertex[3] = {0};
      double decayToTranslated[3] = {0};
      double translationAngle = 0., magDecayToVertex = 0., magDecayToTranslated = 0.; // Use double in case translationAngle is so tiny then gives nan results
      for (int dim = 0; dim < 3; dim++) {
        decayToVertex[dim] = OnAxisVertex[dim]-decaypos[dim];
        decayToTranslated[dim] = OffAxisVertex[dim]-decaypos[dim];

        translationAngle += (decayToVertex[dim])*(decayToTranslated[dim]);
        magDecayToVertex += pow(decayToVertex[dim], 2);
        magDecayToTranslated += pow(decayToTranslated[dim], 2);
      }
      magDecayToVertex = sqrt(magDecayToVertex);
      magDecayToTranslated = sqrt(magDecayToTranslated);
      translationAngle /= (magDecayToVertex*magDecayToTranslated);
      translationAngle = acos(translationAngle);

      // Calculate rotation axis
      // Cross-product
      float translationAxis[3] = {0};
      translationAxis[0] = decayToVertex[1]*decayToTranslated[2] - decayToVertex[2]*decayToTranslated[1];
      translationAxis[1] = decayToVertex[2]*decayToTranslated[0] - decayToVertex[0]*decayToTranslated[2];
      translationAxis[2] = decayToVertex[0]*decayToTranslated[1] - decayToVertex[1]*decayToTranslated[0];
      float magTranslationAxis = 0.;
      for (int dim = 0; dim < 3; dim++) magTranslationAxis += pow(translationAxis[dim], 2);
      magTranslationAxis = sqrt(magTranslationAxis);
      // Get rotation axis n_{hat}
      if(magTranslationAxis!=0)  {for (int dim = 0; dim < 3; dim++) translationAxis[dim] /= magTranslationAxis;}
      else{for (int dim = 0; dim < 3; dim++) {translationAxis[dim] = 0.;translationAngle =0.;}}
      Eigen::Affine3f rTranslation_NDtoND(Eigen::Affine3f(Eigen::AngleAxisf(translationAngle, Eigen::Vector3f(translationAxis[0], translationAxis[1], translationAxis[2]))));

      // Combine: we only have translation when moving ND onaxis to offaxis
      rThrow_NDtoND = rTranslation_NDtoND;

      if(verbosity)
      {
        std::cout << "translationAngle: " << translationAngle << std::endl;
        std::cout << "translationAxis: " << translationAxis[0] << ", " << translationAxis[1] << ", " << translationAxis[2] << "\n";
      }
    }

    // Put everything together in single transform and store.
    // transforms_NDtoND.emplace_back(tThrow_NDtoND * tBack_NDtoND * rThrow_NDtoND * tThere_NDtoND);
    transforms_NDtoND.emplace_back(tBack_NDtoND * rThrow_NDtoND * tThere_NDtoND );

    /*
    I want to apply a rotation to an event and then move it to a different place.
    First I move the event vertex to the origin of the coordinate system with tThere.
    Then I apply the rotation, rThrow. Since I moved the event vertex to the origin, I know the rotation will not move the vertex, as intended.
    Then, move the vertex back to its original position, tBack.
    And finally, move it to the new position with tThrow
    */

    return transforms_NDtoND;
    // returen value: 1D vector
}

// Set Sim_mu_end_vertex
void geoEff::setMuEndV(float x, float y, float z){
  OffAxisMuEndV_BF.resize(3);
  OffAxisMuEndV_BF.at(0)=x;
  OffAxisMuEndV_BF.at(1)=y;
  OffAxisMuEndV_BF.at(2)=z;
}

// Get Sim_mu_end_vertex after rotations
float geoEff::getOffAxisMuEndV(int dim){

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > VectorCoordinate(OffAxisMuEndV_BF.data(),3,OffAxisMuEndV_BF.size()/3,Eigen::OuterStride<>(3));
  // Get the rotated vector coordinate, a 3*1 matrix
  Eigen::Matrix3Xf RotMuEndV_AF = getTransforms_NDtoND()[0] * VectorCoordinate;
  // Return the results for (x,y,z)<->dim=(0,1,2)
  return RotMuEndV_AF(dim, 0);
}

// Set ND_Sim_mu_hadronic_hit
void geoEff::setHadronHitV(float x, float y, float z){
  OffAxisHadronHitV_BF.resize(3);
  OffAxisHadronHitV_BF.at(0)=x;
  OffAxisHadronHitV_BF.at(1)=y;
  OffAxisHadronHitV_BF.at(2)=z;
}

// Get Sim_mu_end_vertex after rotations
float geoEff::getOffAxisHadronHitV(int dim){

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > VectorCoordinate(OffAxisHadronHitV_BF.data(),3,OffAxisHadronHitV_BF.size()/3,Eigen::OuterStride<>(3));
  // Get the rotated vector coordinate, a 3*1 matrix
  Eigen::Matrix3Xf OffAxisHadronHitV_AF = getTransforms_NDtoND()[0] * VectorCoordinate;
  // Return the results for (x,y,z)<->dim=(0,1,2)
  return OffAxisHadronHitV_AF(dim, 0);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Momentum vector space Eigen transformation
std::vector< Eigen::Transform<float,3,Eigen::Affine> > geoEff::getTransforms_NDtoND_P(){

  std::vector< Eigen::Transform<float,3,Eigen::Affine> > transforms_NDtoND_P;

    // Rotation
    Eigen::Affine3f rThrow_NDtoND_P;
    {
      // Calculate rotation due to translation
      // Calculate rotation angle
      double decayToVertex[3] = {0};
      double decayToTranslated[3] = {0};
      double translationAngle = 0, magDecayToVertex = 0, magDecayToTranslated = 0;
      for (int dim = 0; dim < 3; dim++) {
        decayToVertex[dim] = OnAxisVertex[dim]-decaypos[dim];
        decayToTranslated[dim] = OffAxisVertex[dim]-decaypos[dim];

        translationAngle += (decayToVertex[dim])*(decayToTranslated[dim]);
        magDecayToVertex += pow(decayToVertex[dim], 2);
        magDecayToTranslated += pow(decayToTranslated[dim], 2);
      }
      magDecayToVertex = sqrt(magDecayToVertex);
      magDecayToTranslated = sqrt(magDecayToTranslated);
      translationAngle /= (magDecayToVertex*magDecayToTranslated);
      translationAngle = acos(translationAngle);

      // Calculate rotation axis
      // Cross-product
      float translationAxis[3] = {0};
      translationAxis[0] = decayToVertex[1]*decayToTranslated[2] - decayToVertex[2]*decayToTranslated[1];
      translationAxis[1] = decayToVertex[2]*decayToTranslated[0] - decayToVertex[0]*decayToTranslated[2];
      translationAxis[2] = decayToVertex[0]*decayToTranslated[1] - decayToVertex[1]*decayToTranslated[0];
      float magTranslationAxis = 0.;
      for (int dim = 0; dim < 3; dim++) magTranslationAxis += pow(translationAxis[dim], 2);
      magTranslationAxis = sqrt(magTranslationAxis);
      // Get rotation axis n_{hat}
      if(magTranslationAxis!=0)  {for (int dim = 0; dim < 3; dim++) translationAxis[dim] /= magTranslationAxis;}
      else{for (int dim = 0; dim < 3; dim++) {translationAxis[dim] = 0.;translationAngle =0.;}}


      Eigen::Affine3f rTranslation_NDtoND_P(Eigen::Affine3f(Eigen::AngleAxisf(translationAngle, Eigen::Vector3f(translationAxis[0], translationAxis[1], translationAxis[2]))));

      rThrow_NDtoND_P = rTranslation_NDtoND_P;
    }

    // Put everything together in single transform and store.
    // Momentum space is different to the position space, it is not relevant to the position of a vector, so we just need rThrow.
    transforms_NDtoND_P.emplace_back(rThrow_NDtoND_P);

    return transforms_NDtoND_P;
    // returen value: 1D vector
}


// Set ND_Sim_mu_start_p
void geoEff::setMuStartP(float x, float y, float z){
  OffAxisMuStartP_BF.resize(3);
  OffAxisMuStartP_BF.at(0)=x;
  OffAxisMuStartP_BF.at(1)=y;
  OffAxisMuStartP_BF.at(2)=z;
}

// Get Sim_mu_end_vertex after rotations
float geoEff::getOffAxisMuStartP(int dim){

  // Set the Eigen map
  Eigen::Map<Eigen::Matrix3Xf,0,Eigen::OuterStride<> > VectorCoordinate(OffAxisMuStartP_BF.data(),3,OffAxisMuStartP_BF.size()/3,Eigen::OuterStride<>(3));
  // Get the rotated vector coordinate
  Eigen::Matrix3Xf RotMuStartP_AF = getTransforms_NDtoND_P()[0] * VectorCoordinate;

  return RotMuStartP_AF(dim, 0);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Vector doesn't change
double geoEff::RemainUnchanged(double element)
{
  return element;
}
// Momentum calculations
float geoEff::getTotalMomentum(double momentum[3])
{
  float TotalP = sqrt(pow(momentum[0],2)+pow(momentum[1],2)+pow(momentum[2],2));
  return TotalP;
}
double geoEff::getDistance(double v1[3],double v2[3])
{
  double distance = sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2)+pow(v1[2]-v2[2],2));
  return distance;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Earth curvature rotations
// Local y-z axes in FD and ND are rotated due to Earth curvature, x direction is not change
// FD event coordinates, if unchanged, would represent different event in ND coordinate sys.
// Apply an active transformation matrix R_x(theta): rotate each point counterclockwise by theta around x-axis in ND
// theta is 2*|beamLineRotation|
// Transform FD relative coordinate, x coordinate unchanged
//
//              [ 1          0             0
// R_x(theta) =   0      cos(theta)   -sin(theta)
//                0      sin(theta)    cos(theta) ]

double geoEff::getEarthCurvature(double v[3], double BeamAngle, int dim)
{
  double Vector_af[3];
  Vector_af[0]=v[0];
  Vector_af[1]=cos( 2*abs(BeamAngle) )*v[1] - sin( 2*abs(BeamAngle) )*v[2];
  Vector_af[2]=sin( 2*abs(BeamAngle) )*v[1] + cos( 2*abs(BeamAngle) )*v[2];
  return Vector_af[dim];
}

// Reimplemented from getEarthCurvature but with rotation reversed
//                      [ 1          0             0
// R_x(theta)_inverse =   0      cos(theta)    sin(theta)
//                        0     -sin(theta)    cos(theta) ]
Eigen::Matrix3Xf geoEff::getn2fEarthCurvatureCorr(Eigen::Matrix3Xf EdepsposMatrix, double BeamAngle)
{
  // Calculate rotation due to earth curvature
  Eigen::Affine3f rECC;
  rECC = Eigen::Affine3f(Eigen::AngleAxisf(-2*abs(BeamAngle), Eigen::Vector3f(1,0,0))); // X-axis is rotation axis
  Eigen::Matrix3Xf ECCposEdeps = rECC * EdepsposMatrix;
  return ECCposEdeps;
}

Eigen::Matrix3Xf geoEff::move2ndorigin(Eigen::Matrix3Xf randndhitSegPosMatrix)
{
  // Move vertex to ND det coordinate system origin
  Eigen::Affine3f t2ndorig(Eigen::Translation3f(Eigen::Vector3f(-ndrandvertex[0], -ndrandvertex[1], -ndrandvertex[2])));

  Eigen::Matrix3Xf vtxNDoriginEdepspos = t2ndorig * randndhitSegPosMatrix;

  return vtxNDoriginEdepspos;
}

struct throwcombo geoEff::moveBack2ndVertex(Eigen::Matrix3Xf randndhitSegPosMatrix, double BeamAngle)
{
  struct throwcombo ndeccthrowcombo;

  // Figure out how many multiples of 64 bits needed to store output
  int n_longs = N_THROWS_NDECC / 64;
  if (N_THROWS_NDECC % 64) n_longs++;

  std::vector< Eigen::Matrix3Xf > transformedEdepssNDECC;

  // Pass/fail for each set of vetoSize and vetoEnergy
  std::vector< std::vector< std::vector< uint64_t > > > NDECCContainment4RandomThrow(vetoSize.size(), std::vector< std::vector< uint64_t > >(vetoEnergy.size(), std::vector < uint64_t >(n_longs, 0)));

  // Calculate rotation due to randomly thrown angle at FD
  // this is a rotation w.r.t. the beam direction at FD {0., sin(|beamLineRotation|), cos(|beamLineRotation|)]
  Eigen::Affine3f rFDbeamThrow;
  rFDbeamThrow = Eigen::Affine3f(Eigen::AngleAxisf(ndeccrotations[0], Eigen::Vector3f(0, sin(abs(BeamAngle)), cos(abs(BeamAngle)))));

  // Move vertex from (0,0,0) back to the random thrown position in ND
  Eigen::Affine3f tBack2ndVertex(Eigen::Translation3f(Eigen::Vector3f(ndecctranslations[0][0], ndecctranslations[1][0], ndecctranslations[2][0])));

  Eigen::Matrix3Xf vtxNDEdepspos = tBack2ndVertex * rFDbeamThrow * randndhitSegPosMatrix;
  transformedEdepssNDECC.emplace_back(vtxNDEdepspos);

  // Here I am cheating as we only have N_THROWS_NDECC = 1 each time, t == 0 is always true
  // So no loop over throws
  for (unsigned int i = 0; i < vetoSize.size(); i++){
    for (unsigned int j = 0; j < vetoEnergy.size(); j++){
      // Check containment and set bit
      if (isContainedInND(vtxNDEdepspos, hitSegEdeps, vetoSize[i], vetoEnergy[j]))
      {
        // Here I am cheating as we only have N_THROWS_NDECC = 1 each time, t == 0 is always true
        NDECCContainment4RandomThrow[i][j][0/64] |= ((uint64_t)1)<<(0%64);
	    }
    }
  }

  ndeccthrowcombo.thrownEdepspos = transformedEdepssNDECC;
  ndeccthrowcombo.containresult = NDECCContainment4RandomThrow;

  return ndeccthrowcombo;
}

// Put events back to beam center
double geoEff::getTranslations(double v_bf[3], double vtx_bf[3], double vtx_af[3], int dim)
{
  double Vector_af[3];
  for(int i=0; i<3; i++)
  {
    Vector_af[i] = v_bf[i] + (vtx_af[i]-vtx_bf[i]);
  }
  return Vector_af[dim];
}

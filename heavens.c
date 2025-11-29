
double G = .01; //gravitational constant

//toapproximate as hyperbolic to prevent weird behavior around e=1
#define HYPERBOLIC_THRESHOLD (0.9999999999999999999)

typedef double vec2[2];

typedef struct Orbit {
    double semiMajorAxis;
    double eccentricity;
    double meanAnomaly;

    //angle deviation
    double periapsisCos; 
    double periapsisSin;

    int orbitBody; //orbiting body (as array index)
    double mu; //standard gravitational parameter
    double period;
    double meanMotion; //average speed
} Orbit;

void printOrbit(Orbit* orbit) {
    printf("a: %f, e: %f, M: %f\nmu: %f P: %f\nmeanMotion: %f body: %d\n", orbit->semiMajorAxis, orbit->eccentricity, orbit->meanAnomaly, orbit->mu, orbit->period, orbit->meanMotion, orbit->orbitBody);
}

void rotate(double* ret, double coss, double sinn) {
    double a = ret[0] * coss - ret[1] * sinn;
    ret[1] = ret[0] * sinn + ret[1] * coss;
    ret[0] = a;
}

double meanAnomaly(Orbit* orbit, double time) {
    return orbit->eccentricity < HYPERBOLIC_THRESHOLD ? M_PI * 2 / orbit->period * time : orbit->meanMotion * time;
}

double eccentricAnomaly(Orbit* orbit) {
    double approx;
    //approximate with newton's method
    if (orbit->eccentricity < HYPERBOLIC_THRESHOLD) {
        approx = orbit->meanAnomaly;
        for (int i = 0; i < 4; i++) {
            approx -= (approx - orbit->eccentricity * sin(approx) - orbit->meanAnomaly) / (1 - orbit->eccentricity * cos(approx));
        }
    } else {
        approx = (orbit->meanAnomaly > 0) ? (.975 * log(orbit->meanAnomaly)) : (-0.975 * log(-orbit->meanAnomaly));
        for (int i = 0; i < 7; i++) {
            approx -= (orbit->eccentricity * sinh(approx) - approx - orbit->meanAnomaly) / (orbit->eccentricity * cosh(approx) - 1);
        }
    }
    return approx;
}

double eccentricAnomaly2(Orbit* orbit, double meanAnomaly) {
    double approx;
    //approximate with newton's method
    if (orbit->eccentricity < HYPERBOLIC_THRESHOLD) {
        approx = meanAnomaly;
        for (int i = 0; i < 4; i++) {
            approx -= (approx - orbit->eccentricity * sin(approx) - meanAnomaly) / (1 - orbit->eccentricity * cos(approx));
        }
    } else {
        approx = (orbit->meanAnomaly > 0) ? (.975 * log(orbit->meanAnomaly)) : (-0.975 * log(-orbit->meanAnomaly));
        for (int i = 0; i < 7; i++) {
            approx -= (orbit->eccentricity * sinh(approx) - approx - meanAnomaly) / (orbit->eccentricity * cosh(approx) - 1);
        }
    }
    return approx;
}

//convert to cartesian coordinates
void orbitPosition(double* ret, Orbit* orbit, double eccentricAnomaly) {
    if (orbit->eccentricity < HYPERBOLIC_THRESHOLD) {
        ret[0] = orbit->semiMajorAxis * (cos(eccentricAnomaly) - orbit->eccentricity);
        ret[1] = orbit->semiMajorAxis * sqrt(1 - orbit->eccentricity * orbit->eccentricity) * sin(eccentricAnomaly);
    } else {
        ret[0] = -orbit->semiMajorAxis * (orbit->eccentricity - cosh(eccentricAnomaly));
        ret[1] = -orbit->semiMajorAxis * sqrt(orbit->eccentricity * orbit->eccentricity - 1) * sinh(eccentricAnomaly);
    }
    rotate(ret, orbit->periapsisCos, orbit->periapsisSin);
}

//time passed, not time total.
void updateMeanAnomaly(Orbit* orbit, double time) {
    orbit->meanAnomaly = (orbit->eccentricity < HYPERBOLIC_THRESHOLD) ? fmod(orbit->meanAnomaly + orbit->meanMotion * time, M_PI * 2) : (orbit->meanAnomaly + orbit->meanMotion * time);
} 

Orbit createOrbit(double sma, double ecc, int orbitBody, double orbitBodyMass, double anglePeriapsis) {
    Orbit a;
    a.semiMajorAxis = sma;
    a.eccentricity = ecc;
    a.meanAnomaly = 0;
    a.periapsisSin = sin(anglePeriapsis);
    a.periapsisCos = cos(anglePeriapsis);

    a.orbitBody = orbitBody;
    a.mu = G * orbitBodyMass;
    if (ecc < HYPERBOLIC_THRESHOLD) {
        a.period = M_PI * 2 * sqrt(sma*sma*sma / a.mu);
        a.meanMotion = M_PI * 2 / a.period;
    } else {
        a.period = 3600;
        a.meanMotion = sqrt(a.mu / abs(sma*sma*sma));
    }

    return a;
}

typedef struct Planet {
    Orbit orbit;
    float color[4];
    float size; //radius
    double mass;
    vec2 pos;
    vec2 vel; //for objects in gravity well. 
    double soiRadius;
    double hillSoiRadius;
} Planet;

double soiRadius(double semiMajorAxis, double mass, double orbitingMass) {
    return semiMajorAxis * pow(mass/orbitingMass, .4);
}

double hillSoiRadius(double semiMajorAxis, double ecc, double mass, double orbitingMass) {
    return semiMajorAxis * (1 - ecc) * pow(mass/(3.0*orbitingMass),1.0/3.0);
}

Planet createPlanet(double sma, double ecc, double periapsisAngle, float size, float mass, int orbitBody, double orbitBodyMass, float colorRed, float colorGreen, float colorBlue, float colorAlpha) {
    Planet a;
    a.color[0] = colorRed;
    a.color[1] = colorGreen;
    a.color[2] = colorBlue;
    a.color[3] = colorAlpha;
    a.vel[0] = 0;
    a.vel[1] = 0;
    a.size = size;
    a.mass = mass;
    if (orbitBody > -1) {
        a.orbit = createOrbit(sma, ecc, orbitBody, orbitBodyMass, periapsisAngle);
        orbitPosition(a.pos, &a.orbit, eccentricAnomaly(&a.orbit));
    } else {
        a.pos[0]=0;
        a.pos[1]=0;
    }
    //laplace
    if (orbitBody == -1 | orbitBodyMass < .00000000001)
        orbitBodyMass = .00001;
    a.soiRadius = soiRadius(sma, mass, orbitBodyMass);
    a.hillSoiRadius = hillSoiRadius(sma, ecc, mass, orbitBodyMass);
    return a;
}

//get array of points for opengl or something start point is at start and end. misleading name, also draws hyperbolic trajectories.
void ellipseCoords(float* ret, Orbit* orbit, int size) {
    double timeTick = orbit->period / (double)(size - 1);
    
    for (int i = 0; i < size; i++) {
        double pos[2];
        double m = meanAnomaly(orbit, (i - size / 2) * timeTick);
        double E = eccentricAnomaly2(orbit, m);
        orbitPosition(pos, orbit, E);
        
        ret[i*2] = (float)pos[0];
        ret[i*2+1] = (float)pos[1];
        //printf("%f %f %f %f\n",m, E, ret[i*2], ret[i*2+1]);
    }
}

double radialDistance(Orbit* orbit) {
    return orbit->semiMajorAxis * (1 - orbit->eccentricity * cos(eccentricAnomaly(orbit)));
}

double orbitalVelocity(Orbit* orbit, double radialDistance) {
    return sqrt(orbit->mu * (2/radialDistance - 1/orbit->semiMajorAxis));
}

double orbitalVelocity2(Orbit* orbit) {
    double radialDist = radialDistance(orbit);
    return sqrt(orbit->mu * (2/radialDist - 1/orbit->semiMajorAxis));
}

//the z value number in the cross product (x and y are 0 in 2d)
double crossProduct1(double a1, double a2, double b1, double b2) {
    return a1*b2 - b1*a2;
}

//returns 2d vector (x and y are )
void crossProduct2(double* ret, double a1, double a2, double h) {
    ret[0] = a2*h;
    ret[1] = -a1*h;
}

//must be given velocities relevant to orbiting body
Orbit createOrbitFromVelocity(double x, double y, double xv, double yv, int orbitBody, double orbitBodyMass) {
    double v = sqrt(xv*xv + yv*yv);
    double r = sqrt(x*x + y*y);
    Orbit a;
    a.orbitBody = orbitBody;
    a.mu = G * orbitBodyMass;
    double specificOrbitalEnergy = v*v/2 - a.mu / r;
    a.semiMajorAxis = -a.mu / (2 * specificOrbitalEnergy);
    double h = crossProduct1(x, y, xv, yv);

    double ecc[2];
    crossProduct2(ecc, xv, yv, h);
    ecc[0] = ecc[0] / a.mu - (x / r);
    ecc[1] = ecc[1] / a.mu - (y / r);

    a.eccentricity = sqrt(ecc[0]*ecc[0]+ecc[1]*ecc[1]);
    double trueAnomaly = acos((ecc[0] * x + ecc[1] * y) / (a.eccentricity * r));
    if (x*xv + y*yv < 0)
        trueAnomaly += M_PI;


    //vector towards periapsis
    a.periapsisCos = ecc[0] / a.eccentricity;
    a.periapsisSin = ecc[1] / a.eccentricity;
    if (a.eccentricity < HYPERBOLIC_THRESHOLD) {
        double E = 2 * atan(sqrt((1 - a.eccentricity) / (1 + a.eccentricity)) * tan(trueAnomaly/2));
        a.meanAnomaly = E - a.eccentricity * sin(E);
        a.period = M_PI * 2 * sqrt(a.semiMajorAxis*a.semiMajorAxis*a.semiMajorAxis / a.mu);
        a.meanMotion = M_PI * 2 / a.period;
    } else {
        //double E = 2 * atanh(sqrt((1 + a.eccentricity) / (a.eccentricity - 1)) * tan(trueAnomaly/2));
        double p = r * (1 + a.eccentricity * cos(trueAnomaly));
        double E = asinh(sqrt(a.mu / abs(p)) * tan(trueAnomaly) * sqrt(a.eccentricity - 1));
        a.meanAnomaly = a.eccentricity * sinh(E) - E;
        a.period = 3600;
        a.meanMotion = sqrt(a.mu / abs(a.semiMajorAxis*a.semiMajorAxis*a.semiMajorAxis));
    }

    return a;
}

typedef struct Ship {
    Orbit orbit;
    vec2 pos;
    vec2 vel;
    double orient;
    double mass;
    bool orbitChanged;
    bool warping;
} Ship;

Ship createShip(double x, double y, double xv, double yv, double orient, double mass, int orbitBody, double orbitBodyMass) {
    Ship a;
    a.pos[0] = x;
    a.pos[1] = y;
    a.vel[0] = xv;
    a.vel[1] = yv;
    a.orient = orient;
    a.orbit = createOrbitFromVelocity(x, y, xv, yv, orbitBody, orbitBodyMass);
    a.orbitChanged = true;
    a.mass = mass;
    a.warping = false;

    return a;
}

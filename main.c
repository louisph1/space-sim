#ifdef WIN32
void* __chk_fail=0; //janky fix
#endif

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "heavens.c"

#define CIRCLE_SIZE (128)

uint32_t VBO;
uint32_t VAO;
uint32_t EBO;
unsigned int shaderProgram;

uint32_t windowWidth=600;
uint32_t windowHeight=600;

//mouse data
float mx=0.001f;
float my=0.001f;
bool mdown = 0;
bool m2down = 0;

bool projectionChanged = true;

const char* getGLErrorStr(GLenum err)
{
    switch (err)
    {
    case GL_NO_ERROR:          return "No error";
    case GL_INVALID_ENUM:      return "Invalid enum";
    case GL_INVALID_VALUE:     return "Invalid value";
    case GL_INVALID_OPERATION: return "Invalid operation";
    //case GL_STACK_OVERFLOW:    return "Stack overflow";
    //case GL_STACK_UNDERFLOW:   return "Stack underflow";
    case GL_OUT_OF_MEMORY:     return "Out of memory";
    default:                   return "Unknown error";
    }
}

float randf() {
    return (float)rand() / (float)RAND_MAX;
}

int mod (int a, int b) {
    int c = a % b;
    return c < 0 ? c + b : c;
}

void mouse_callback(GLFWwindow* window, double dmposx, double dmposy) {
    mx = (float)(dmposx) - (float)(windowWidth/2);
    my = -1.0f*((float)(dmposy) - (float)(windowHeight/2));
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        mdown = 1;
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        mdown = 0;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        m2down = 1;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
        m2down = 0;
    }
}

void window_size_callback(GLFWwindow* window, int width, int height) {
    glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
    glViewport(0,0,windowWidth,windowHeight);
    projectionChanged = true;
}

int main() {
    srand (time(NULL));
    
    const char* vertexShaderSrc = "#version 330 core\n"
    "uniform mat4 projection;"
    "uniform vec2 camPos;"
    "uniform float size;"
    "layout (location=0) in vec2 triPos;"
    "void main(){"
    "vec2 npos=triPos*size - camPos;"
    "gl_Position=projection * vec4(npos,0,1);}\0";

    const char* fragmentShaderSrc = "#version 330 core\n"
    "out vec4 outColor;"
    "uniform vec4 color;"
    "void main(){outColor=vec4(color);}\0";

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(windowWidth,windowHeight,"test",NULL,NULL);
    if(window==NULL){
        printf("no window :(");
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        printf("Failed to initialize GLAD\n");
        return -1;
    }   
    
    //capture mouse
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED); 
    glfwSetCursorPosCallback(window, mouse_callback);  
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    
    glViewport(0,0,windowWidth,windowHeight);
    
    
    //shaders
    int success;
    char infoLog[512];
    uint32_t vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader,1,&vertexShaderSrc,NULL);
    glCompileShader(vertexShader);

    glGetShaderiv(vertexShader,GL_COMPILE_STATUS,&success);
    if(!success){
        glGetShaderInfoLog(vertexShader,512,NULL,infoLog);
        printf(infoLog);
        return 1;
    }

    uint32_t fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader,1,&fragmentShaderSrc,NULL);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader,GL_COMPILE_STATUS,&success);
    if(!success){
        glGetShaderInfoLog(fragmentShader,512,NULL,infoLog);
        printf(infoLog);
        return 1;
    }

    shaderProgram=glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram,GL_LINK_STATUS,&success);
        
    if(!success){
        glGetProgramInfoLog(shaderProgram,512,NULL,infoLog);
        printf(infoLog);
        return 1;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    glUseProgram(shaderProgram);
    
    glGenVertexArrays(1,&VAO);
    glGenBuffers(1,&VBO);
    glGenBuffers(1,&EBO);

    glBindVertexArray(VAO);   
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*512, NULL, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 128, NULL, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,2*sizeof(float),(void*)0);
    glEnableVertexAttribArray(0);








    //init graphics
    const int camPosLocation = glGetUniformLocation(shaderProgram, "camPos");
    const int colorLocation = glGetUniformLocation(shaderProgram, "color");
    const int projectionLocation = glGetUniformLocation(shaderProgram, "projection");
    const int sizeLocation = glGetUniformLocation(shaderProgram, "size");
    
    float aspect = (float)windowWidth/(float)windowHeight;
    //mat4 projection;
    //glm_ortho(-1.0f*windowWidth/2,windowWidth/2,-1.0f*windowHeight/2,windowHeight/2,-1,1, projection);
    
    float zoom = 10;
    float invZoom = 1 / zoom;

    float projection[] = {
        invZoom/windowWidth, 0.000000, 0.000000, 0.000000,
        0.000000, invZoom/windowHeight, 0.000000, 0.000000,
        0.000000, 0.000000, -1.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 1.000000
    };
    
    glUniformMatrix4fv(projectionLocation, 1, 0, projection);

    //transparency
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glfwSwapInterval(1);

    float circleDraw[CIRCLE_SIZE * 2];
    circleDraw[0] = 0.0f;
    circleDraw[1] = 0.0f;

    for (int i=1;i<CIRCLE_SIZE;i++) {
        circleDraw[i*2] = cos((float)i * M_PI * 2.0f / ((float)CIRCLE_SIZE-2.0f));
        circleDraw[i*2+1] = sin((float)i * M_PI * 2.0f / ((float)CIRCLE_SIZE-2.0f));
    }

    float shipVerts[] = {
        1.0f,0.0f,
        cos(M_PI*5/6), -sin(M_PI*5/6),
        cos(M_PI*5/6), sin(M_PI*5/6)
    };

    float camx = 0.0f;
    float camy = 0.0f;
    float realCamx = 0.0f;
    float realCamy = 0.0f;
    
    //init game

    int numPlanets = 8;
    int trackedPlanet = numPlanets;

    Planet planets[numPlanets];
    //Planet createPlanet(double sma, double ecc, double angle, float size, float mass, int orbitBody, double orbitBodyMass, float colorRed, float colorGreen, float colorBlue, float colorAlpha) {
    planets[0] = createPlanet(0, 0, 0, 50, 1000000, -1, 0, 1, 1, 0, 1); //sun
    planets[1] = createPlanet(200, 0.22, M_PI_2, 4, 4, 0, 1000000, .7, .35, .15, 1); //mercury
    planets[2] = createPlanet(1000, 0.1, 0, 10, 1000, 0, 1000000, 0, 0, 1, 1); //earth
    planets[3] = createPlanet(36, 0.03, 0, 2, 10, 2, 100, .8, .8, .8, 1); //moon
    planets[4] = createPlanet(1200, 0.05, 0, 8, 800, 0, 1000000, .8, .3, .1, 1); //mars
    planets[5] = createPlanet(1000, 0.9, 1, .5, 1, 0, 1000000, .5, 5, .5, 1); //asteroid
    planets[6] = createPlanet(700, .9, 0, 8, 180, 0, 1000000, .9, .3, .9, 1); //planet y (test)
    planets[7] = createPlanet(500, .9, 0, 8, 180, 0, 1000000, .9, .3, .9, 1); //planet x (hyperbolic)
    planets[6].pos[0] = 500;
    planets[6].pos[1] = 0;
    planets[6].orbit = createOrbitFromVelocity(500,0,0,4,0,1000000);
    planets[6].soiRadius = soiRadius(planets[6].orbit.semiMajorAxis, planets[6].mass, 1000000);
    planets[6].hillSoiRadius = hillSoiRadius(planets[6].orbit.semiMajorAxis, planets[6].orbit.eccentricity, planets[6].mass, 1000000);
    planets[7].pos[0] = 500;
    planets[7].pos[1] = 100;
    planets[7].orbit = createOrbitFromVelocity(500,100,0,10,0,1000000);
    printOrbit(&planets[7].orbit);

    int orbitDrawingSize = 256;
    float planetOrbitPoints[numPlanets][orbitDrawingSize*2];

    for (int i = 1; i < numPlanets; i++) {
        ellipseCoords(&planetOrbitPoints[i][0], &planets[i].orbit, orbitDrawingSize);
    }

    int numShips = 1;
    Ship ships[numShips];

    ships[0] = createShip(400, 1, 1, 6, 0.5, .0000001, 0, 1000000);

    float shipOrbitPoints[numShips][orbitDrawingSize*2];
    for (int i = 0; i < numShips; i++) {
        ellipseCoords(&shipOrbitPoints[i][0], &ships[i].orbit, orbitDrawingSize);
    }
    








    bool zpressed = false;
    bool xpressed = false;
    bool cpressed = false;
    bool vpressed = false;
    bool orbitView = true;
    bool soiView = true;

    printf("Initialization Successful\n");

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
            realCamy += zoom * 15.0f;
        }
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
            realCamx -= zoom * 15.0f;
        }
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
            realCamy -= zoom * 15.0f;
        }
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
            realCamx += zoom * 15.0f;
        }
        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            zoom *= 1.1f;
            invZoom = 1 / zoom;
            projectionChanged = true;
        }
        if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
            zoom /= 1.1f;
            invZoom = 1 / zoom;
            projectionChanged = true;
        }
        /*if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
            rotate(&ships[0].orbit.periapsisCos, cos(.1),sin(.1));
            ellipseCoords(&shipOrbitPoints[0][0], &ships[0].orbit, orbitDrawingSize);
        }*/


        
        {
            int a = glfwGetKey(window, GLFW_KEY_C);
            if (a == GLFW_PRESS) {
                if (!cpressed) {
                    soiView = !soiView;
                }
                cpressed = true;
            } else if (a == GLFW_RELEASE) {
                cpressed = false;
            }
        }
        {
            int a = glfwGetKey(window, GLFW_KEY_V);
            if (a == GLFW_PRESS) {
                if (!vpressed) {
                    orbitView = !orbitView;
                }
                vpressed = true;
            } else if (a == GLFW_RELEASE) {
                vpressed = false;
            }
        }
        {
            int a = glfwGetKey(window, GLFW_KEY_Z);
            if (a == GLFW_PRESS) {
                if (!zpressed) {
                    trackedPlanet = mod(trackedPlanet-1, numPlanets+1);
                    realCamx = 0;
                    realCamy = 0;
                }
                zpressed = true;
            } else if (a == GLFW_RELEASE) {
                zpressed = false;
            }
        }
        {
            int a = glfwGetKey(window, GLFW_KEY_X);
            if (a == GLFW_PRESS) {
                if (!xpressed) {
                    trackedPlanet = mod(trackedPlanet+1, numPlanets+1);
                    realCamx = 0;
                    realCamy = 0;
                }
                xpressed = true;
            } else if (a == GLFW_RELEASE) {
                xpressed = false;
            }
        }
        







        //game logic
        double timeTick = 1.0/ 30.0;

        //planetary motion
        for (int i = 1; i < numPlanets; i++) {
            updateMeanAnomaly(&planets[i].orbit, timeTick);
            double E = eccentricAnomaly(&planets[i].orbit);
            double oldx = planets[i].pos[0];
            double oldy = planets[i].pos[1];
            orbitPosition(planets[i].pos, &planets[i].orbit, E);
            planets[i].pos[0] += planets[planets[i].orbit.orbitBody].pos[0];
            planets[i].pos[1] += planets[planets[i].orbit.orbitBody].pos[1];
            planets[i].vel[0] = (planets[i].pos[0] - oldx) / timeTick;
            planets[i].vel[1] = (planets[i].pos[1] - oldy) / timeTick;
        }

        //player controls
        double mouseDist = sqrt(mx*mx+my*my);
        double t1=mx / mouseDist;
        double t2=my / mouseDist;
             
        ships[0].orient = atan2(t2,t1);
        
        if (mdown) {
            ships[0].vel[0] += t1 * 0.1 * timeTick;
            ships[0].vel[1] += t2 * 0.1 * timeTick;
            ships[0].orbitChanged = true;
        }

        if (m2down) {
            ships[0].vel[0] += t1 * 5.0 * timeTick;
            ships[0].vel[1] += t2 * 5.0 * timeTick;
            ships[0].orbitChanged = true;
        }
        
        //ship logic
        for (int i = 0; i < numShips; i++) {
            //printf("%f %f\n ", sqrt(ships[i].vel[0]*ships[i].vel[0]+ships[i].vel[1]*ships[i].vel[1]), orbitalVelocity(&planets[1].orbit, radialDistance(&planets[1].orbit)));
            int p = ships[i].orbit.orbitBody;

            int precision = 128;
            for (int j = 0; j < precision; j++) {
                ships[i].pos[0] += ships[i].vel[0] * timeTick / precision / 2;
                ships[i].pos[1] += ships[i].vel[1] * timeTick / precision / 2;
                
                
                int q = p;
                while (1) {
                    double dx = ships[i].pos[0] - planets[q].pos[0];
                    double dy = ships[i].pos[1] - planets[q].pos[1];
                    double r2 = dx*dx+dy*dy;
                    double r = sqrt(r2);
                    double force = G * planets[q].mass / r2;

                    ships[i].vel[0] -= dx / r * force * timeTick / precision;
                    ships[i].vel[1] -= dy / r * force * timeTick / precision;

                    if (q == 0)
                        break;
                    q = planets[q].orbit.orbitBody;
                    
                }

                ships[i].pos[0] += ships[i].vel[0] * timeTick / precision / 2;
                ships[i].pos[1] += ships[i].vel[1] * timeTick / precision / 2;
            }

            double dx = ships[i].pos[0] - planets[p].pos[0];
            double dy = ships[i].pos[1] - planets[p].pos[1];
                
            if (p > 0 && dx*dx+dy*dy > planets[p].soiRadius*planets[p].soiRadius) {
                int n = planets[p].orbit.orbitBody;
                ships[i].orbit = createOrbitFromVelocity(dx, dy, ships[i].vel[0]-planets[0].vel[0], ships[i].vel[1]-planets[0].vel[1], n, planets[n].mass);
            }

            for (int j = p + 1; j < numPlanets; j++) {
                double a = ships[i].pos[0] - planets[j].pos[0];
                double b = ships[i].pos[1] - planets[j].pos[1];
                if (a*a+b*b < planets[j].soiRadius*planets[j].soiRadius) {
                    ships[i].orbit = createOrbitFromVelocity(dx, dy, ships[i].vel[0]-planets[j].vel[0], ships[i].vel[1]-planets[j].vel[1], j, planets[j].mass);
                }
            }
        }





        //render
        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);

        if (projectionChanged) {
            projection[0] = 1.0f * invZoom / windowWidth;
            projection[5] = 1.0f * invZoom / windowHeight;
            glUniformMatrix4fv(projectionLocation, 1, 0, projection);
            projectionChanged = false;
        }
        
        //planet tracking
        if (trackedPlanet == numPlanets) {
            camx = realCamx + ships[0].pos[0];
            camy = realCamy + ships[0].pos[1];
        } else {
            camx = realCamx + planets[trackedPlanet].pos[0];
            camy = realCamy + planets[trackedPlanet].pos[1];
        }
        
        //ship orbital drawings
        for (int i = 0; i < numShips; i++) {
            if (ships[i].orbitChanged) {
                int p = ships[i].orbit.orbitBody;
                double dx = ships[i].pos[0] - planets[p].pos[0];
                double dy = ships[i].pos[1] - planets[p].pos[1];
                ships[i].orbit = createOrbitFromVelocity(dx,dy,ships[i].vel[0]-planets[p].vel[0],ships[i].vel[1]-planets[p].vel[1], p, planets[p].mass);
                ellipseCoords(&shipOrbitPoints[i][0], &ships[i].orbit, orbitDrawingSize);
                ships[i].orbitChanged = false;
            }
        }
        
        glClearColor(0.01f, 0.01f, 0.01f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        //views
        if (soiView) {
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(circleDraw), &circleDraw);
            for (int i = 0; i < numPlanets; i++) {
                glUniform1f(sizeLocation, planets[i].soiRadius);
                glUniform2f(camPosLocation, camx-planets[i].pos[0], camy-planets[i].pos[1]);
                glUniform4f(colorLocation, 1,1,1,.4);
                glDrawArrays(GL_TRIANGLE_FAN,0,CIRCLE_SIZE);

                //hill soi (more orbitally stable, for looks)
                glUniform1f(sizeLocation, planets[i].hillSoiRadius);
                glUniform4f(colorLocation, 1,1,1,.1);
                glDrawArrays(GL_TRIANGLE_FAN,0,CIRCLE_SIZE);
            }
        }

        if (orbitView) {
            //draw orbit lines
            glUniform1f(sizeLocation, 1);
            glUniform4f(colorLocation, 1,1,1,1);
            for (int i = 1; i < numPlanets; i++) {
                glUniform2f(camPosLocation, camx-planets[planets[i].orbit.orbitBody].pos[0], camy-planets[planets[i].orbit.orbitBody].pos[1]);
                glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*orbitDrawingSize*2, &planetOrbitPoints[i][0]);
                glDrawArrays(GL_LINE_STRIP,0,orbitDrawingSize);
            }
            for (int i = 0; i < numShips; i++) {
                glUniform2f(camPosLocation, camx-planets[ships[i].orbit.orbitBody].pos[0], camy-planets[ships[i].orbit.orbitBody].pos[1]);
                glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*orbitDrawingSize*2, &shipOrbitPoints[i][0]);
                glDrawArrays(GL_LINE_STRIP,0,orbitDrawingSize);
            }
        }

        //draw planets
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(circleDraw), &circleDraw);
        for (int i = 0; i < numPlanets; i++) {
            glUniform1f(sizeLocation, planets[i].size);
            glUniform2f(camPosLocation, camx-planets[i].pos[0], camy-planets[i].pos[1]);
            glUniform4f(colorLocation, planets[i].color[0],planets[i].color[1],planets[i].color[2],planets[i].color[3]);
            glDrawArrays(GL_TRIANGLE_FAN,0,CIRCLE_SIZE);
        }

        //draw ships
        glUniform4f(colorLocation, 1,1,1,1);
        glUniform1f(sizeLocation, zoom * 15);
        float adjShipVerts[6];
        for (int i = 0; i < numShips; i++) {
            float sina = sin(ships[i].orient);
            float cosa = cos(ships[i].orient);
            for (int i = 0; i < 6; i+=2) {
                adjShipVerts[i] = shipVerts[i]*cosa - shipVerts[i+1]*sina;
                adjShipVerts[i+1] = shipVerts[i]*sina + shipVerts[i+1]*cosa;
            }
            glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(adjShipVerts), &adjShipVerts);
            glUniform2f(camPosLocation, camx-ships[i].pos[0], camy-ships[i].pos[1]);
            glDrawArrays(GL_TRIANGLES,0,3);
        }

        
        glfwSwapBuffers(window);
        
        //error test
        GLenum e = glGetError();
        if (e) {
            printf(getGLErrorStr(e));
            fflush(stdout);
        }
        
        usleep(1000000/60); //60fps
    }
    
    printf("bye\n");
    
    return 0;
}
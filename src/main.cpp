#define GLM_ENABLE_EXPERIMENTAL


#include <GL3/gl3.h>
#include <GL3/gl3w.h>
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/ext/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale
#include <glm/ext/matrix_clip_space.hpp> // glm::perspective
#include <glm/ext/scalar_constants.hpp> // glm::pi
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/intersect.hpp>

#include <GLFW/glfw3.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <random>

#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"





#define BUFFER_OFFSET(a) ((void*)(a))
#define MAX_MATERIAL_COUNT 100

#ifndef GLFW_TRUE
#define GLFW_TRUE 1
#endif

#ifndef GLFW_FALSE
#define GLFW_FALSE 0
#endif

#define PI 3.14159265358979323846


typedef struct {
	GLenum       type;
	const char* filename;
	GLuint       shader;
} ShaderInfo;


struct Triangle {
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 v2;
	glm::vec3 Norm[3];
	glm::vec3 face_normal;
	unsigned char Color[3];
};

struct BoundingBox {
	glm::vec3 min;
	glm::vec3 max;
};

struct Ray {
    glm::vec3 origin;
    glm::vec3 direction;
    float energy;
    int bounceCount;
};

struct Intersection {
    bool hit;
    float distance;
    glm::vec3 point;
    glm::vec3 normal;
    int triangleIndex;
};

struct Reflection {
    float delay;
    float amplitude;
    glm::vec3 direction;
    int bounceCount; 
};

// WAV file header structure
#pragma pack(push, 1)
struct WAVHeader {
    // RIFF chunk
    char riffId[4] = {'R', 'I', 'F', 'F'};
    uint32_t riffSize;
    char waveId[4] = {'W', 'A', 'V', 'E'};
    
    // fmt chunk
    char fmtId[4] = {'f', 'm', 't', ' '};
    uint32_t fmtSize = 16;
    uint16_t audioFormat = 1; // PCM
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
    
    // data chunk
    char dataId[4] = {'d', 'a', 't', 'a'};
    uint32_t dataSize;
};
#pragma pack(pop)

int g_IRSampleRate = 44100;
float g_IRDuration = 3.0f; // in seconds
bool g_NormalizeIR = true;

//enum VAO_IDs { Triangles, NumVAOs };
enum Buffer_IDs { VtxBuffer = 0, NormBuffer = 1, ColBuffer = 2, NumBuffers = 3 };
enum Attrib_IDs { vPosition = 0, vNormalVertex = 1, vColor = 2};

enum RenderMode { RENDER_SOLID, RENDER_WIREFRAME, RENDER_POINTS };
RenderMode g_RenderMode = RENDER_SOLID;

enum WindingMode { CW, CCW };
WindingMode g_WindingMode = CCW;

enum CameraMode { LOOKAT_CAM, FREE_CAM };
CameraMode g_CameraMode = LOOKAT_CAM;

enum ShadingMode { NO_SHADING, GOURAUD, PHONG };
ShadingMode g_ShadingMode = PHONG;

GLuint  VAO;
GLuint  Buffers[NumBuffers];
GLuint  NumVertices;

GLuint g_ColorBufferTexture = 0;
GLuint g_FramebufferQuadVAO = 0;
GLuint g_FramebufferQuadVBO = 0;
GLuint g_FramebufferProgram = 0;

float g_AspectRatio = 1.0f; //actual value will be calculated in FramebufferSizeCallback
int g_WindowWidth = 800;
int g_WindowHeight = 600;


// Rendering stuff
float constexpr max_FOV = 90.0f;
float constexpr min_FOV = 5.0f;
float constexpr default_FOV = 60.0f;
// os valores padr�o aqui s�o modificados em tempo de execu��o
glm::vec3 g_CameraPos(0.0f, 0.0f, 0.0f); // coordenadas cartesianas da c�mera no sistema de coordenadas global
//coordenadas esf�ricas da c�mera
float g_CameraTheta = 0.0f; // �ngulo no plano ZX em rela��o ao eixo Z 
float g_CameraPhi = 0.0f;   // �ngulo em rela��o ao eixo Y 
float g_CameraDistance = 6.0f; // Dist�ncia da c�mera para o ponto lookat
float g_FOV = default_FOV; // degrees 
float g_nearplane = 0.1f;  // Dist�ncia do "near plane"
float g_farplane = 3000.0f; // Dist�ncia do "far plane"
glm::vec3 g_lookat_point(0.0f, 0.0f, 0.0f); // Ponto para onde a c�mera est� olhando
BoundingBox g_bbox; // Bounding box do modelo
glm::vec3 g_color(1.0f, 1.0f, 1.0f); // Cor usada para renderizar o modelo
glm::vec3 g_LightDir = glm::normalize(glm::vec3(1.0f, -2.0f, 1.0f));


// Audio ray-tracing stuff
float g_RoomScale = 1.0f;
int g_NumRays = 1000;
int g_MaxBounces = 3;
float g_ListenerRadius = 0.1f;
float g_DirectTime = 0.0f;
glm::vec3 g_ListenerPosition(0.0f,0.0f,0.0f);
glm::vec3 g_SourcePosition(0.0f,0.0f,1.0f);
std::vector<Reflection> g_Reflections;
std::vector<Triangle> g_SceneTriangles; // Yeah, another structure to store triangle data




GLint g_view_matrix_uniform, g_model_matrix_uniform, g_proj_matrix_uniform;


bool g_LeftMousePressed = false;
bool g_RightMousePressed = false;
bool g_MiddleMousePressed = false;
double g_LastMousePosX, g_LastMousePosY;


static const GLchar* ReadShader(const char* filename);
GLuint LoadShaders(ShaderInfo* shaders);
int readModel(char* filename, float** Vert, float** Vert_Normal);
BoundingBox getBoundingBox(float* Vert, int NumTris);
glm::vec3 getBoundingBoxCenter(const BoundingBox& bbox);
float calculateCameraDistance(const BoundingBox& bbox, float fovY);
glm::vec3 initializeCameraPosition(const glm::vec3& modelCenter, float cameraDistance);
void resetCamera();
void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void CursorPositionCallback(GLFWwindow* window, double xpos, double ypos);
void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);
void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
void FramebufferSizeCallback(GLFWwindow* window, int width, int height);
void ErrorCallback(int error, const char* description);

void traceRayBounces(Ray ray);
void traceAudioRays();
bool rayReachesListener(const glm::vec3& rayPos, const glm::vec3& rayDir, float maxDistance);
glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal) ;
Intersection intersectRayWithScene(const Ray& ray) ;
void buildSceneTriangles(Triangle* tris, int numTris); 

void renderIR();
void writeWAVFile(const std::vector<float>& samples, const std::string& filename, 
                 int sampleRate, int numChannels, int bitsPerSample);
std::string generateFilename();



int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
		return 1;
	}

	char* filename = argv[1];

	glfwInit();

	glfwSetErrorCallback(ErrorCallback);

	GLFWwindow* window = glfwCreateWindow(g_WindowWidth, g_WindowHeight, "CMP143", NULL, NULL);

	glfwMakeContextCurrent(window);
	gl3wInit();


	// backface culling
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	// depth testing
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	// callbacks
	glfwSetMouseButtonCallback(window, MouseButtonCallback);
	glfwSetCursorPosCallback(window, CursorPositionCallback);
	glfwSetScrollCallback(window, ScrollCallback);
	glfwSetKeyCallback(window, KeyCallback);
	glfwSetFramebufferSizeCallback(window, FramebufferSizeCallback);

    FramebufferSizeCallback(window, g_WindowWidth, g_WindowHeight); 

	ShaderInfo  shaders[] =
	{
		{ GL_VERTEX_SHADER, "./triangles.vert" },
		{ GL_FRAGMENT_SHADER, "./triangles.frag" },
		{ GL_NONE, NULL }
	};

	GLuint program = LoadShaders(shaders);


	glUseProgram(program);

	// Initialize ImGui
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 130");



	//////////////////////////////////////////////////////////////////////////

	float* Vert = NULL; 
	float* Vert_Normal = NULL;
	int NumTris;

	// Carregar dados do modelo
	NumTris = readModel(filename, &Vert, &Vert_Normal);
	NumVertices = NumTris * 3;
	std::vector<glm::vec4> Vert_Vec4;
	for (int i = 0; i < NumTris * 3; i++) {
		Vert_Vec4.push_back(glm::vec4(Vert[3 * i], Vert[3 * i + 1], Vert[3 * i + 2], 1.0f));
	}
	std::vector<glm::vec3> Norm_Vec3;
	for (int i = 0; i < NumTris * 3; i++) {
		Norm_Vec3.push_back(glm::vec3(Vert_Normal[3 * i], Vert_Normal[3 * i + 1], Vert_Normal[3 * i + 2]));
	}
	g_bbox = getBoundingBox(Vert, NumTris);
	glm::vec3 modelCenter = getBoundingBoxCenter(g_bbox);
	g_lookat_point = glm::vec3(0.0f,0.0f,0.0f);
	float initialCameraDistance = calculateCameraDistance(g_bbox, g_FOV);
	g_CameraDistance = initialCameraDistance;
	g_CameraPos = initializeCameraPosition(glm::vec3(0.0f), g_CameraDistance);

	glm::mat4 translate_to_center = glm::translate(glm::mat4(1.0f), -modelCenter);


	// Store triangles that will be used for ray tracing
	Triangle* triangles = new Triangle[NumTris];

	// temporary hack, proper way should be to translate the model to the origin
	g_ListenerPosition = glm::vec3(modelCenter.x, modelCenter.y, modelCenter.z - (abs(g_bbox.min.z - g_bbox.max.z)*0.25) );
	g_SourcePosition = glm::vec3(modelCenter.x, modelCenter.y, modelCenter.z + (abs(g_bbox.min.z - g_bbox.max.z)*0.25) );

	//print some useful information
	std::cout << "Model center: (" 
			  << modelCenter.x << ", " 
			  << modelCenter.y << ", " 
			  << modelCenter.z << ")" << std::endl;
	std::cout << "Bounding Box Min: (" 
			  << g_bbox.min.x << ", " << g_bbox.min.y << ", " << g_bbox.min.z << ")\n";
	std::cout << "Bounding Box Max: (" 
			  << g_bbox.max.x << ", " << g_bbox.max.y << ", " << g_bbox.max.z << ")\n";
	std::cout << "Source position: (" 
			  << g_SourcePosition.x << ", " 
			  << g_SourcePosition.y << ", " 
			  << g_SourcePosition.z << ")" << std::endl;
	std::cout << "Listener position: (" 
			  << g_ListenerPosition.x << ", " 
			  << g_ListenerPosition.y << ", " 
			  << g_ListenerPosition.z << ")" << std::endl;

			  
	//calculate listener and source positions
	float directDistance = glm::length(g_SourcePosition - g_ListenerPosition);
	float speedOfSound = 343.0f; // m/s
	float directTime = directDistance / speedOfSound;
	g_DirectTime = directTime;

	//more useful information
	std::cout << "Direct source-listener distance: " << directDistance << " meters" << std::endl;
	std::cout << "Direct sound travel time: " << directTime << " seconds" << std::endl;


	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	glCreateBuffers(NumBuffers, Buffers);
	glBindBuffer(GL_ARRAY_BUFFER, Buffers[VtxBuffer]);
	glBufferStorage(GL_ARRAY_BUFFER, Vert_Vec4.size() * 4 * sizeof(GL_FLOAT), Vert_Vec4.data(), GL_DYNAMIC_STORAGE_BIT);
	glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
	glEnableVertexAttribArray(vPosition);
	glBindBuffer(GL_ARRAY_BUFFER, Buffers[NormBuffer]);
	glBufferStorage(GL_ARRAY_BUFFER, NumTris * 9 * sizeof(GL_FLOAT), Vert_Normal, GL_DYNAMIC_STORAGE_BIT);
	glVertexAttribPointer(vNormalVertex, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0));
	glEnableVertexAttribArray(vNormalVertex);
	GLint model_matrix_uniform = glGetUniformLocation(program, "model");
	GLint proj_matrix_uniform = glGetUniformLocation(program, "projection");
	GLint view_matrix_uniform = glGetUniformLocation(program, "view");
	GLint color_uniform = glGetUniformLocation(program, "in_color");
	GLint light_dir_uniform = glGetUniformLocation(program, "lightDir");
	GLint shading_mode_uniform = glGetUniformLocation(program, "shadingMode"); // Add this line
	glm::mat4 view_matrix;
	glm::mat4 projection_matrix;
	glm::mat4 model_matrix;
	glm::mat4 MVP_matrix;

	while (!glfwWindowShouldClose(window))
	{
		static const float black[] = { 0.0f, 0.0f, 0.0f, 0.0f };
		glClearBufferfv(GL_COLOR, 0, black);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		switch (g_RenderMode) {
		case RENDER_SOLID:
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			break;
		case RENDER_WIREFRAME:
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			break;
		case RENDER_POINTS:
			glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
			break;
		default:
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		switch (g_WindingMode) {
		case CW:
			glFrontFace(GL_CW);
			break;
		case CCW:
			glFrontFace(GL_CCW);
			break;
		default:
			glFrontFace(GL_CCW);
		}
				
		// rendering is the same as in the programming assignments
		if (g_CameraMode == LOOKAT_CAM) {
			g_CameraPos.y = g_CameraDistance * sin(g_CameraPhi) + modelCenter.y;
			g_CameraPos.z = (g_CameraDistance * cos(g_CameraPhi) * cos(g_CameraTheta)) + modelCenter.z;
			g_CameraPos.x = (g_CameraDistance * cos(g_CameraPhi) * sin(g_CameraTheta)) + modelCenter.x;
			glm::vec4 camera_position_c = glm::vec4(g_CameraPos.x, g_CameraPos.y, g_CameraPos.z, 1.0f); 
			glm::vec4 camera_lookat_l = glm::vec4(g_lookat_point, 1.0f); 
			glm::vec4 camera_up_vector = glm::vec4(0.0f, 1.0f, 0.0f, 0.0f);
			view_matrix = glm::lookAt(glm::vec3(camera_position_c), glm::vec3(camera_lookat_l), glm::vec3(camera_up_vector));
		}
		else if (g_CameraMode == FREE_CAM) {
			glm::vec3 cam_forward = -glm::normalize(glm::vec3(
				cos(g_CameraPhi) * sin(g_CameraTheta),
				sin(g_CameraPhi),
				cos(g_CameraPhi) * cos(g_CameraTheta)
			));
			glm::vec3 cam_right = glm::normalize(glm::vec3(
				sin(g_CameraPhi - glm::pi<float>() * 0.5f),
				0,
				cos(g_CameraPhi - glm::pi<float>() * 0.5f)
			));
			glm::vec3 cam_up = glm::cross(cam_forward, cam_right);

			view_matrix = glm::lookAt(g_CameraPos, g_CameraPos + cam_forward, cam_up);
		}
		projection_matrix = glm::perspective(glm::radians(g_FOV), g_AspectRatio, g_nearplane, g_farplane);
		model_matrix = translate_to_center;
		glUseProgram(program);
		glUniformMatrix4fv(view_matrix_uniform, 1, GL_FALSE, glm::value_ptr(view_matrix));
		glUniformMatrix4fv(proj_matrix_uniform, 1, GL_FALSE, glm::value_ptr(projection_matrix));
		glUniformMatrix4fv(model_matrix_uniform, 1, GL_FALSE, glm::value_ptr(model_matrix));
		glUniform3fv(color_uniform, 1, glm::value_ptr(g_color));
		glUniform3fv(light_dir_uniform, 1, glm::value_ptr(g_LightDir));
		glUniform1i(shading_mode_uniform, static_cast<int>(g_ShadingMode)); // Add this line
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, Buffers[VtxBuffer]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, Vert_Vec4.size() * 4 * sizeof(GL_FLOAT), Vert_Vec4.data());
		glDrawArrays(GL_TRIANGLES, 0, NumVertices);
	
		// Imgui controls
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::Begin("Controls");

		// ray tracing controls
		ImGui::Text("Audio Ray Tracing");
		if (ImGui::Button("Trace Audio Rays")) {
			traceAudioRays();
		}
		if (ImGui::Button("Render IR")) {
			renderIR();
		}
		ImGui::SliderInt("Number of Rays", &g_NumRays, 10, 100000);
		ImGui::SliderInt("Max Bounces", &g_MaxBounces, 1, 20);
		ImGui::SliderFloat("Room Scale", &g_RoomScale, 0.1f, 50.0f);
		ImGui::SliderFloat("Listener Radius", &g_ListenerRadius, 0.01f, 2.0f);

		//camera and rendering controls
		ImGui::Text("Camera Controls");
		if (ImGui::Button("Reset Camera")) {
			resetCamera();
		}
        ImGui::Text("Camera Mode");  
        ImGui::RadioButton("Look-At Camera", (int*)&g_CameraMode, LOOKAT_CAM); ImGui::SameLine();
        ImGui::RadioButton("Free Camera", (int*)&g_CameraMode, FREE_CAM);
		ImGui::SliderFloat("Distance", &g_CameraDistance, 1.0f, 2500);
		ImGui::Text("Winding Mode:");
		ImGui::RadioButton("CCW", (int*)&g_WindingMode, CCW); ImGui::SameLine();
		ImGui::RadioButton("CW", (int*)&g_WindingMode, CW);
		ImGui::Text("Render Mode:");
		ImGui::RadioButton("Solid", (int*)&g_RenderMode, RENDER_SOLID); ImGui::SameLine();
		ImGui::RadioButton("Wireframe", (int*)&g_RenderMode, RENDER_WIREFRAME); ImGui::SameLine();

		ImGui::End();
		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);

	glfwTerminate();
}



static const GLchar* ReadShader(const char* filename)
{
	FILE* infile = fopen(filename, "rb");

	if (!infile) {
		std::cerr << "Unable to open file '" << filename << "'" << std::endl;
		return NULL;
	}

	fseek(infile, 0, SEEK_END);
	int len = ftell(infile);
	fseek(infile, 0, SEEK_SET);

	GLchar* source = new GLchar[len + 1];

	fread(source, 1, len, infile);
	fclose(infile);

	source[len] = 0;

	return const_cast<const GLchar*>(source);
}

GLuint LoadShaders(ShaderInfo* shaders)
{
	if (shaders == NULL) { return 0; }

	GLuint program = glCreateProgram();

	ShaderInfo* entry = shaders;
	while (entry->type != GL_NONE) {
		GLuint shader = glCreateShader(entry->type);

		entry->shader = shader;

		const GLchar* source = ReadShader(entry->filename);
		if (source == NULL) {
			for (entry = shaders; entry->type != GL_NONE; ++entry) {
				glDeleteShader(entry->shader);
				entry->shader = 0;
			}

			return 0;
		}

		glShaderSource(shader, 1, &source, NULL);
		delete[] source;

		glCompileShader(shader);

		GLint compiled;
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
		if (!compiled) {
			GLsizei len;
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);

			GLchar* log = new GLchar[len + 1];
			glGetShaderInfoLog(shader, len, &len, log);
			std::cerr << "Shader compilation failed: " << log << std::endl;
			delete[] log;

			return 0;
		}

		glAttachShader(program, shader);

		++entry;
	}

	glLinkProgram(program);

	GLint linked;
	glGetProgramiv(program, GL_LINK_STATUS, &linked);
	if (!linked) {
		GLsizei len;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);

		GLchar* log = new GLchar[len + 1];
		glGetProgramInfoLog(program, len, &len, log);
		std::cerr << "Shader linking failed: " << log << std::endl;
		delete[] log;

		for (entry = shaders; entry->type != GL_NONE; ++entry) {
			glDeleteShader(entry->shader);
			entry->shader = 0;
		}

		return 0;
	}

	return program;
}

int readModel(char* filename, float** Vert, float** Vert_Normal) {
	glm::vec3 ambient[MAX_MATERIAL_COUNT], diffuse[MAX_MATERIAL_COUNT], specular[MAX_MATERIAL_COUNT];
	float shine[MAX_MATERIAL_COUNT];

	int NumTris, material_count, color_index[3], i;
	char ch;

	Triangle* Tris;


	FILE* fp = fopen(filename, "r");
	if (fp == NULL) { printf("ERROR: unable to open TriObj [%s]!\n", filename); return -1; }

	fscanf(fp, "%c", &ch);
	while (ch != '\n') // skip the first line � object�s name 
		fscanf(fp, "%c", &ch);

	fscanf(fp, "# triangles = %d\n", &NumTris);     // read # of triangles 
	fscanf(fp, "Material count = %d\n", &material_count);    // read material count

	for (i = 0; i < material_count; i++) {
		fscanf(fp, "ambient color %f %f %f\n", &(ambient[i].x), &(ambient[i].y), &(ambient[i].z));
		fscanf(fp, "diffuse color %f %f %f\n", &(diffuse[i].x), &(diffuse[i].y), &(diffuse[i].z));
		fscanf(fp, "specular color %f %f %f\n", &(specular[i].x), &(specular[i].y), &(specular[i].z));
		fscanf(fp, "material shine %f\n", &(shine[i]));
	}

	fscanf(fp, "%c", &ch);
	while (ch != '\n')      // skip documentation line 
		fscanf(fp, "%c", &ch);
	// 
	//  allocate triangles for tri model 
	// 
	printf("Reading in %s (%d triangles). . .\n", filename, NumTris);
	Tris = new Triangle[NumTris];
	// 
	for (i = 0; i < NumTris; i++)     // read triangles 
	{
		fscanf(fp, "v0 %f %f %f %f %f %f %d\n",
			&(Tris[i].v0.x), &(Tris[i].v0.y), &(Tris[i].v0.z),
			&(Tris[i].Norm[0].x), &(Tris[i].Norm[0].y), &(Tris[i].Norm[0].z),
			&(color_index[0]));
		fscanf(fp, "v1 %f %f %f %f %f %f %d\n",
			&(Tris[i].v1.x), &(Tris[i].v1.y), &(Tris[i].v1.z),
			&(Tris[i].Norm[1].x), &(Tris[i].Norm[1].y), &(Tris[i].Norm[1].z),
			&(color_index[1]));
		fscanf(fp, "v2 %f %f %f %f %f %f %d\n",
			&(Tris[i].v2.x), &(Tris[i].v2.y), &(Tris[i].v2.z),
			&(Tris[i].Norm[2].x), &(Tris[i].Norm[2].y), &(Tris[i].Norm[2].z),
			&(color_index[2]));
		fscanf(fp, "face normal %f %f %f\n", &(Tris[i].face_normal.x), &(Tris[i].face_normal.y),
			&(Tris[i].face_normal.z));
		// 
		Tris[i].Color[0] = (unsigned char)(int)(255 * (diffuse[color_index[0]].x));
		Tris[i].Color[1] = (unsigned char)(int)(255 * (diffuse[color_index[0]].y));
		Tris[i].Color[2] = (unsigned char)(int)(255 * (diffuse[color_index[0]].z));
	}
	fclose(fp);

	buildSceneTriangles(Tris, NumTris);

	// 
	//     For use in the vertex buffer objects in your application, pack the vertex and normal data 
	//           into vectors 
	// 
	*Vert = new float[9 * NumTris];
	*Vert_Normal = new float[9 * NumTris];

	for (i = 0; i < NumTris; i++) {
		//    vertex coordinates  
		(*Vert)[9 * i] = Tris[i].v0.x;
		(*Vert)[9 * i + 1] = Tris[i].v0.y;
		(*Vert)[9 * i + 2] = Tris[i].v0.z;
		(*Vert)[9 * i + 3] = Tris[i].v1.x;
		(*Vert)[9 * i + 4] = Tris[i].v1.y;
		(*Vert)[9 * i + 5] = Tris[i].v1.z;
		(*Vert)[9 * i + 6] = Tris[i].v2.x;
		(*Vert)[9 * i + 7] = Tris[i].v2.y;
		(*Vert)[9 * i + 8] = Tris[i].v2.z;

		//    vertex normal coordinates  

		(*Vert_Normal)[9 * i] = Tris[i].Norm[0].x;
		(*Vert_Normal)[9 * i + 1] = Tris[i].Norm[0].y;
		(*Vert_Normal)[9 * i + 2] = Tris[i].Norm[0].z;
		(*Vert_Normal)[9 * i + 3] = Tris[i].Norm[1].x;
		(*Vert_Normal)[9 * i + 4] = Tris[i].Norm[1].y;
		(*Vert_Normal)[9 * i + 5] = Tris[i].Norm[1].z;
		(*Vert_Normal)[9 * i + 6] = Tris[i].Norm[2].x;
		(*Vert_Normal)[9 * i + 7] = Tris[i].Norm[2].y;
		(*Vert_Normal)[9 * i + 8] = Tris[i].Norm[2].z;
	}
	return NumTris;
}

// Calcula a "bounding box" do modelo
BoundingBox getBoundingBox(float* Vert, int NumTris) {
	BoundingBox bbox;
	bbox.min = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	bbox.max = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	for (int i = 0; i < NumTris * 3; i++) {
		glm::vec3 vertex = glm::vec3(Vert[3 * i], Vert[3 * i + 1], Vert[3 * i + 2]);
		bbox.min = glm::min(bbox.min, vertex);
		bbox.max = glm::max(bbox.max, vertex);
	}

	return bbox;
}

glm::vec3 getBoundingBoxCenter(const BoundingBox& bbox) {
	return (bbox.min + bbox.max) * 0.5f;

}

// Calcula a dist�ncia necess�ria para que o modelo caiba na tela
float calculateCameraDistance(const BoundingBox& bbox, float fovY) {
	glm::vec3 bboxSize = bbox.max - bbox.min;
	float bboxDiagonal = glm::length(bboxSize);
	float distance = (bboxDiagonal / 2.0f) / tan(glm::radians(fovY) / 2.0f);
	return distance;
}

// inicializa a posic�o da c�mera de forma que o modelo caiba na tela e esteja posicionado no centro
glm::vec3 initializeCameraPosition(const glm::vec3& modelCenter, float cameraDistance) {
	glm::vec3 cameraPosition = modelCenter + glm::vec3(0.0f, 0.0f, cameraDistance);
	return cameraPosition;
}

// reseta a posi��o e orienta��o da c�mera aos valores padr�o
void resetCamera() {
	g_FOV = default_FOV;
	g_CameraDistance = calculateCameraDistance(g_bbox, g_FOV);
	g_CameraPos = initializeCameraPosition(g_lookat_point, g_CameraDistance);
	g_CameraTheta = 0.0f;
	g_CameraPhi = 0.0f;
}

void MouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		if (action == GLFW_PRESS) {
			g_LeftMousePressed = true;
			glfwGetCursorPos(window, &g_LastMousePosX, &g_LastMousePosY);
		}
		else if (action == GLFW_RELEASE) {
			g_LeftMousePressed = false;
		}
	}
	else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (action == GLFW_PRESS) {
			g_RightMousePressed = true;
			glfwGetCursorPos(window, &g_LastMousePosX, &g_LastMousePosY);
		}
		else if (action == GLFW_RELEASE) {
			g_RightMousePressed = false;
		}
	}
	else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
		if (action == GLFW_PRESS) {
			g_MiddleMousePressed = true;
			glfwGetCursorPos(window, &g_LastMousePosX, &g_LastMousePosY);
		}
		else if (action == GLFW_RELEASE) {
			g_MiddleMousePressed = false;
		}
	}
}

void CursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
	if (g_LeftMousePressed) {
		double dx = xpos - g_LastMousePosX;
		double dy = ypos - g_LastMousePosY;
		static float sensitivity = 0.005f;

		if (g_CameraMode == FREE_CAM) {
			sensitivity = -0.005f;
		}
		else {
			sensitivity = 0.005f;
		}

		g_CameraTheta += static_cast<float>(dx) * sensitivity;
		g_CameraPhi -= static_cast<float>(dy) * sensitivity;

		if (g_CameraPhi > glm::pi<float>() * 0.5f) {
			g_CameraPhi = glm::pi<float>() * 0.5f;
		}
		if (g_CameraPhi < -glm::pi<float>() * 0.5f) {
			g_CameraPhi = -glm::pi<float>() * 0.5f;
		}
		

		g_LastMousePosX = xpos;
		g_LastMousePosY = ypos;

		//printf("Yaw: %f, Pitch: %f\n", g_CameraYaw, g_CameraPitch);
	}
	if (g_RightMousePressed) { 
		;
	}
}

void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
	g_CameraDistance -= 10.0f * yoffset;
	const float verysmallnumber = std::numeric_limits<float>::epsilon();
	if (g_CameraDistance < verysmallnumber)
		g_CameraDistance = verysmallnumber;
}

void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	static const float moveSpeed = 10.0f;
	if (action == GLFW_PRESS || action == GLFW_REPEAT) {
		if(g_CameraMode == FREE_CAM) {
			glm::vec3 cam_forward = -glm::normalize(glm::vec3(
				cos(g_CameraPhi) * sin(g_CameraTheta),
				sin(g_CameraPhi),
				cos(g_CameraPhi) * cos(g_CameraTheta)
			));
			glm::vec3 cam_right = glm::normalize(glm::vec3(
				sin(g_CameraPhi - glm::pi<float>() * 0.5f),
				0,
				cos(g_CameraPhi - glm::pi<float>() * 0.5f)
			));
			glm::vec3 cam_up = glm::cross(cam_forward, cam_right);

			if (key == GLFW_KEY_UP) {
				g_CameraPos += cam_forward * moveSpeed;
			}
			if (key == GLFW_KEY_DOWN) {
				g_CameraPos -= cam_forward * moveSpeed;
			}
			if (key == GLFW_KEY_LEFT) {
				g_CameraPos -= cam_right * moveSpeed;
			}
			if (key == GLFW_KEY_RIGHT) {
				g_CameraPos += cam_right * moveSpeed;
			}
			if (key == GLFW_KEY_PAGE_DOWN) {
				g_CameraPos -= cam_up * moveSpeed;
			}
			if (key == GLFW_KEY_PAGE_UP) {
				g_CameraPos += cam_up * moveSpeed;
			}
		}
	}
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
	else if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		resetCamera();
	}
	else if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
		g_color = glm::vec3(1.0f, 0.0f, 0.0f); // red
	}
	else if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
		g_color = glm::vec3(0.0f, 1.0f, 0.0f); // green
	}
	else if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
		g_color = glm::vec3(0.0f, 0.0f, 1.0f); // blue
	}
	else if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
		g_color = glm::vec3(1.0f, 1.0f, 0.0f); // yellow
	}
	else if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
		g_color = glm::vec3(1.0f, 0.0f, 1.0f); // magenta
	}
	else if (key == GLFW_KEY_6 && action == GLFW_PRESS) {
		g_color = glm::vec3(0.0f, 1.0f, 1.0f); // cyan
	}
	else if (key == GLFW_KEY_7 && action == GLFW_PRESS) {
		g_color = glm::vec3(1.0f, 1.0f, 1.0f); // white
	}
	else if (key == GLFW_KEY_LEFT_BRACKET && action == GLFW_PRESS) {
		g_nearplane -= 10.0f;
		if (g_nearplane < 0.1f) g_nearplane = 0.1f;
		std::cout << "Near plane changed: " << g_nearplane << std::endl;
	}
	else if (key == GLFW_KEY_RIGHT_BRACKET && action == GLFW_PRESS) {
		g_nearplane += 10.0f;
		std::cout << "Near plane changed: " << g_nearplane << std::endl;
	}
	else if (key == GLFW_KEY_MINUS && action == GLFW_PRESS) {
		g_farplane -= 100.0f;
		if (g_farplane < 10.0f) g_farplane = 10.0f;
		std::cout << "Far plane changed: " << g_farplane << std::endl;
	}
	else if (key == GLFW_KEY_EQUAL && action == GLFW_PRESS) {
		g_farplane += 100.0f;
		std::cout << "Far plane changed: " << g_farplane << std::endl;
	}
	else if (key == GLFW_KEY_COMMA && action == GLFW_PRESS) {
		g_FOV -= 3.0f;
		if (g_FOV < min_FOV) g_FOV = min_FOV;
		std::cout << "FOV changed: " << g_FOV << std::endl;
	}
	else if (key == GLFW_KEY_PERIOD && action == GLFW_PRESS) {
	
		g_FOV += 3.0f;
		if (g_FOV > max_FOV) g_FOV = max_FOV;
		std::cout << "FOV changed: " << g_FOV << std::endl;
	}
	else if (key == GLFW_KEY_W && action == GLFW_PRESS) {
		g_RenderMode = RENDER_WIREFRAME;
	}
	else if (key == GLFW_KEY_P && action == GLFW_PRESS) {
		g_RenderMode = RENDER_POINTS;
	}
	else if (key == GLFW_KEY_S && action == GLFW_PRESS) {
		g_RenderMode = RENDER_SOLID;
	}
	else if (key == GLFW_KEY_C && action == GLFW_PRESS) {
		g_WindingMode = (g_WindingMode == CW) ? CCW : CW;
	}
	else if (key == GLFW_KEY_V && action == GLFW_PRESS) {
		g_CameraMode = (g_CameraMode == LOOKAT_CAM) ? FREE_CAM : LOOKAT_CAM;
	}
}

void FramebufferSizeCallback(GLFWwindow* window, int width, int height)
{
	g_WindowWidth = width;
	g_WindowHeight = height;
    glViewport(0, 0, width, height);
    g_AspectRatio = (float)width / height;

}

void ErrorCallback(int error, const char* description)
{
	std::cerr << "Error: GLFW: " << description << std::endl;
}

// store triangles in new structure for ray-tracing
void buildSceneTriangles(Triangle* tris, int numTris) {
    g_SceneTriangles.clear();
    g_SceneTriangles.reserve(numTris);
    for (int i = 0; i < numTris; i++) {
        g_SceneTriangles.push_back(tris[i]);
    }
}

// Ray-triangle intersection using GLM
Intersection intersectRayWithScene(const Ray& ray) {
    Intersection closestHit;
    closestHit.hit = false;
    closestHit.distance = FLT_MAX;
    
    for (int i = 0; i < g_SceneTriangles.size(); i++) {
        const Triangle& tri = g_SceneTriangles[i];
        
        glm::vec2 baryPosition;
		float distance;
		
		bool hit = glm::intersectRayTriangle(
			ray.origin,
			ray.direction,
			tri.v0, tri.v1, tri.v2,
			baryPosition,
			distance
		);
        
        if (hit && distance > 0.001f && distance < closestHit.distance) {
            closestHit.hit = true;
            closestHit.distance = distance;
            closestHit.point = ray.origin + ray.direction * distance;
			
            //closestHit.normal = tri.face_normal;
			glm::vec3 geometricNormal = glm::normalize(tri.face_normal);
			// Ensure the normal points against the ray direction (towards the incoming ray)
			if (glm::dot(ray.direction, geometricNormal) > 0.0f) {
				geometricNormal = -geometricNormal;
			}
			closestHit.normal = geometricNormal;

            closestHit.triangleIndex = i;
        }
    }
    
    return closestHit;
}

// Calculate reflected ray direction
glm::vec3 reflect(const glm::vec3& incident, const glm::vec3& normal) {
    return incident - 2.0f * glm::dot(incident, normal) * normal;
}

bool rayReachesListener(const glm::vec3& rayPos, const glm::vec3& rayDir, float maxDistance) {
    glm::vec3 listenerPos = g_ListenerPosition;
    
    // Calculate vector to listener
    glm::vec3 toListener = listenerPos - rayPos;
    float distanceToListener = glm::length(toListener);
    
    if (distanceToListener > maxDistance) return false;
    
    // Check if ray is not facing away from listener
    float projDist = glm::dot(toListener, rayDir);
    if (projDist <= 0) return false; // Ray pointing away
    
    // Calculate how close ray passes to listener center
    float distSq = distanceToListener * distanceToListener - projDist * projDist;
    if (distSq > g_ListenerRadius * g_ListenerRadius) return false;
    
    // Check for occlusion between reflection point and listener (very inefficient!)
    Ray testRay;
    testRay.origin = rayPos;
    testRay.direction = glm::normalize(toListener);
    testRay.energy = 1.0f;
    testRay.bounceCount = 0;
    Intersection occlusionHit = intersectRayWithScene(testRay);
    if (occlusionHit.hit && occlusionHit.distance < distanceToListener) {
        return false; 
    }
    
    return true;
}

// Main ray tracing function for audio
void traceAudioRays() {
    g_Reflections.clear();
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
    
    glm::vec3 sourcePos = g_SourcePosition;
    
    for (int rayIndex = 0; rayIndex < g_NumRays; rayIndex++) {
        // Generate random ray direction (uniform sphere distribution)
        glm::vec3 direction;
        do {
            direction = glm::vec3(dis(gen), dis(gen), dis(gen));
        } while (glm::length(direction) > 1.0f);
        direction = glm::normalize(direction);
        
        Ray ray;
        ray.origin = sourcePos;
        ray.direction = direction;
        ray.energy = 1.0f;
        ray.bounceCount = 0;
        
        // Trace this ray through multiple bounces
        traceRayBounces(ray);
    }
    
    std::cout << "Found " << g_Reflections.size() << " reflections" << std::endl;
    for (int i = 0; i < (int)g_Reflections.size(); i++) {
        printf("Reflection %d: Delay=%.8fs, Amplitude=%.3f, Bounces=%d\n", 
               i, g_Reflections[i].delay, g_Reflections[i].amplitude, 
               g_Reflections[i].bounceCount);
    }
}


void traceRayBounces(Ray ray) {
    const float speedOfSound = 343.0f; // m/s
    float totalDistance = 0.0f;
    
    
    for (int bounce = 0; bounce < g_MaxBounces; bounce++) {
        // Find intersection with scene
        Intersection hit = intersectRayWithScene(ray);
        
        if (!hit.hit) { // Ray rides off into infinity (hits nothing)
            break; 
        }
        
        // Apply room scale to the distance for audio calculations
        totalDistance += hit.distance * g_RoomScale;
        
		// TODO: implement attenuation based on distance
        // Simple energy attenuation based on amount of bounces
        ray.energy *= 0.8f; 
        
        if (ray.energy < 0.01f) {
            break; // Ray energy too low
        }
        
        // Create reflected ray
        Ray reflectedRay;
        reflectedRay.origin = hit.point + hit.normal * 0.001f; // Offset to avoid self-intersection
        reflectedRay.direction = reflect(ray.direction, hit.normal); //this is the magic line
        reflectedRay.energy = ray.energy;
        reflectedRay.bounceCount = ray.bounceCount + 1;
        
		if (rayReachesListener(reflectedRay.origin, reflectedRay.direction, 1000.0f)) {
			glm::vec3 listenerPos = g_ListenerPosition;
			float distanceToListener = glm::length(listenerPos - reflectedRay.origin) * g_RoomScale;
			float fullDistance = totalDistance + distanceToListener;
			
			Reflection reflection;
			reflection.delay = fullDistance / speedOfSound;
			reflection.amplitude = reflectedRay.energy;
			reflection.direction = reflectedRay.direction;
			reflection.bounceCount = reflectedRay.bounceCount; 
			
			g_Reflections.push_back(reflection);
			break; //Stop tracing if ray reaches listener (listener occludes sound)
		
        }
        
        ray = reflectedRay;
    }
}

// Generate a timestamp-based filename for the WAV file
std::string generateFilename() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << "ir_" << std::put_time(std::localtime(&time), "%Y%m%d_%H%M%S") << ".wav";
    return ss.str();
}

// Render IR from reflections
void renderIR() {
    if (g_Reflections.empty()) {
        std::cout << "No reflections to render. Run ray tracing first." << std::endl;
        return;
    }

    int numSamples = static_cast<int>(g_IRSampleRate * g_IRDuration);
    std::vector<float> irSamples(numSamples, 0.0f);
    
    // Find maximum delay to ensure all reflections fit in our IR duration
    float maxDelay = 0.0f;
    for (const auto& reflection : g_Reflections) {
        maxDelay = std::max(maxDelay, reflection.delay);
    }
    
    // Add the direct time to get the full delay needed
    float maxTotalDelay = maxDelay + g_DirectTime;
    
    if (maxTotalDelay > g_IRDuration) {
        std::cout << "Warning: Reflections with pre-delay exceed IR duration. "
                  << "Max delay: " << maxTotalDelay << "s, IR duration: " 
                  << g_IRDuration << "s" << std::endl;
    }
    
    std::cout << "Rendering IR with " << g_Reflections.size() 
              << " reflections, applying pre-delay of " << g_DirectTime 
              << "s..." << std::endl;
    
	// this could be improved to have sub-sample accuracy, see Wayverb's Early Reflection implementation
    // Place impulses at appropriate time positions with pre-delay
    for (const auto& reflection : g_Reflections) {
        // Add direct time as pre-delay to each reflection
        float totalDelay = reflection.delay + g_DirectTime;

		
        int sampleIndex = static_cast<int>(totalDelay * g_IRSampleRate);
        
        if (sampleIndex < numSamples) {
            // Apply reflection amplitude to the sample
            irSamples[sampleIndex] += reflection.amplitude;
        }
    }
    
    // Normalize if requested
    if (g_NormalizeIR) {
        float maxAmplitude = 0.0f;
        for (const float sample : irSamples) {
            maxAmplitude = std::max(maxAmplitude, std::abs(sample));
        }
        
        if (maxAmplitude > 0.0f) {
            float normFactor = 0.95f / maxAmplitude; // Leave a bit of headroom
            for (float& sample : irSamples) {
                sample *= normFactor;
            }
        }
    }
    
    // Write WAV file
    std::string filename = generateFilename();
    writeWAVFile(irSamples, filename, g_IRSampleRate, 1, 16);
    std::cout << "IR saved to " << filename << std::endl;
}

// Write WAV file from IR samples
void writeWAVFile(const std::vector<float>& samples, const std::string& filename,
                 int sampleRate, int numChannels, int bitsPerSample) {
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    // Calculate sizes for header
    uint32_t dataSize = static_cast<uint32_t>(samples.size() * numChannels * (bitsPerSample / 8));
    uint32_t riffSize = 36 + dataSize;
    
    // Create and configure WAV header
    WAVHeader header;
    header.riffSize = riffSize;
    header.numChannels = numChannels;
    header.sampleRate = sampleRate;
    header.byteRate = sampleRate * numChannels * (bitsPerSample / 8);
    header.blockAlign = numChannels * (bitsPerSample / 8);
    header.bitsPerSample = bitsPerSample;
    header.dataSize = dataSize;
    
    // Write header
    outFile.write(reinterpret_cast<const char*>(&header), sizeof(WAVHeader));
    
    // Write sample data
    if (bitsPerSample == 16) {
        for (float sample : samples) {
            // Convert float to 16-bit integer
            int16_t intSample = static_cast<int16_t>(sample * 32767.0f);
            outFile.write(reinterpret_cast<const char*>(&intSample), 2);
        }
    } else if (bitsPerSample == 24) {
        for (float sample : samples) {
            // Convert float to 24-bit integer (3 bytes)
            int32_t intSample = static_cast<int32_t>(sample * 8388607.0f);
            char bytes[3] = {
                static_cast<char>(intSample & 0xFF),
                static_cast<char>((intSample >> 8) & 0xFF),
                static_cast<char>((intSample >> 16) & 0xFF)
            };
            outFile.write(bytes, 3);
        }
    } else if (bitsPerSample == 32) {
        for (float sample : samples) {
            // Just write the float directly for 32-bit
            outFile.write(reinterpret_cast<const char*>(&sample), 4);
        }
    }
    
    outFile.close();
}
#ifndef CLOSE2GL_H
#define CLOSE2GL_H

#include <stdexcept> 

#include <glm/glm.hpp>
#include <GL3/gl3.h>
#include "matrices.h"


bool isOrthonormal(const glm::vec4& u, const glm::vec4& v, const glm::vec4& n);


namespace Close2GL {
	
	// Classe para representar c�mera
	class Camera {
	public:
		// Construtor
		Camera(glm::vec4 pos, glm::vec4 u, glm::vec4 v, glm::vec4 n,
			float hfov, float vfov, float near_d, float far_d)
			: pos(pos), u(u), v(v), n(n),
			hfov(hfov), vfov(vfov), near_d(near_d), far_d(far_d) {
			if (!isOrthonormal(u, v, n)) {
				throw std::invalid_argument("u, v, and n do not form an orthonormal basis.");
			}
		}
		void set_LookAtDirection(glm::vec3 lookat_point, glm::vec3 up); //Modifica a orienta��o da c�mera para "olhar" a um determinado ponto
		void rotate_u(float angle); //Rotaciona a c�mera em volta do seu pr�prio eixo u 
		void rotate_v(float angle); //Rotaciona a c�mera em volta do seu pr�prio eixo v

		glm::vec4 pos; //Posi��o da c�mera em coordenadas do mundo
		// Vetores que definem o sistema de coordenadas da c�mera
		glm::vec4 u;
		glm::vec4 v;
		glm::vec4 n;
		//Campo de vis�o vertical e horizontal em graus
		float hfov;
		float vfov;
		//Dist�ncias dos planos near e far
		float near_d;
		float far_d;
	};

	//implementa��o das fun��es est� em Close2GL.cpp
	glm::mat4 ViewMatrix(glm::vec4 u, glm::vec4 v, glm::vec4 n, glm::vec4 pos);
	glm::mat4 ViewMatrix_lookAt(glm::vec4 position_c, glm::vec4 lookat_point, glm::vec4 up_vector);
	glm::mat4 ProjectionMatrix_frustum(float l, float r, float b, float t, float n, float f);
	glm::mat4 ProjectionMatrix_perspective(float fovy, float aspect, float near_d, float far_d);
	glm::mat4 ViewportMatrix(int x, int y, int width, int height);
	
	
}

#endif // CLOSE2GL_H

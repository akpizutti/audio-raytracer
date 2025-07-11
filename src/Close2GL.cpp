#include "Close2GL.h"

//Modifica a orienta��o da c�mera para "olhar" a um determinado ponto
void Close2GL::Camera::set_LookAtDirection(glm::vec3 lookat_point, glm::vec3 up) {
	glm::vec4 view_vector = glm::vec4(lookat_point,1.0f) - this->pos;
	this->n = -1.0f * view_vector;
	this->n.w = 0.0f;
	glm::vec4 u = crossproduct(glm::vec4(up,0.0f), this->n);
	this->n = this->n / norm(this->n);
	this->u = u / norm(u);
	this->v = crossproduct(this->n, this->u);

}

//Rotaciona a c�mera em volta do seu pr�prio eixo u 
void Close2GL::Camera::rotate_u(float angle) {
	glm::mat4 rotationMatrix = Matrix_Rotate(angle, this->u);

	this->v = rotationMatrix * this->v;
	this->n = rotationMatrix * this->n;
}

//Rotaciona a c�mera em volta do seu pr�prio eixo v 
void Close2GL::Camera::rotate_v(float angle) {
	glm::mat4 rotationMatrix = Matrix_Rotate(angle, this->v);

	this->u = rotationMatrix * this->u;
	this->n = rotationMatrix * this->n;
}

//Calcula a matriz View, que transforma um ponto para o sistema de coordenadas da c�mera
glm::mat4 Close2GL::ViewMatrix(glm::vec4 u, glm::vec4 v, glm::vec4 n, glm::vec4 pos) {
	pos.w = 0.0f;
	glm::mat4 view_matrix = Matrix(
		u.x, u.y, u.z, -dotproduct(pos,u),
		v.x, v.y, v.z, -dotproduct(pos,v),
		n.x, n.y, n.z, -dotproduct(pos,n),
		0.0f, 0.0f, 0.0f, 1.0f
	);
	return view_matrix;
}

//Calcula a matriz View para uma c�mera do tipo look-at
glm::mat4 Close2GL::ViewMatrix_lookAt(glm::vec4 position_c, glm::vec4 lookat_point, glm::vec4 up_vector) {
	glm::vec4 view_vector = lookat_point - position_c;
	glm::vec4 n = -view_vector;
	glm::vec4 u = crossproduct(up_vector, n);
	n = n / norm(n);
	u = u / norm(u);
	glm::vec4 v = crossproduct(n, u);
	return Close2GL::ViewMatrix(u, v, n, position_c);
}

//Calcula a matriz de proje��o para um view frustum arbitr�rio
glm::mat4 Close2GL::ProjectionMatrix_frustum(float l, float r, float b, float t, float n, float f) {
	glm::mat4 proj_matrix = Matrix(
		(2 * n) / (r - l), 0.0f, (r + l) / (r - l), 0.0f,
		0.0f, (2 * n) / (t - b), (t + b) / (t - b), 0.0f,
		0.0f, 0.0f, -(f + n) / (f - n), -(2 * f * n) / (f - n),
		0.0f, 0.0f, -1.0f, 0.0f
	);
	return proj_matrix;
}

//Calcula a matriz de proje��o para um view frustum sim�trico
glm::mat4 Close2GL::ProjectionMatrix_perspective(float fovy, float aspect, float near_d, float far_d) {
   float top = near_d * tan(fovy / 2.0f);
   float bottom = -top;
   float right = top * aspect;
   float left = -right;

   return Close2GL::ProjectionMatrix_frustum(left, right, bottom, top, near_d, far_d);
}

//Calcula a matriz viewport 
glm::mat4 Close2GL::ViewportMatrix(int x, int y, int width, int height) {
	glm::mat4 viewport_matrix = Matrix(
		width / 2.0f, 0.0f, 0.0f, x + (width / 2.0f),
		0.0f, height / 2.0f, 0.0f, y + (height / 2.0f),
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
	return viewport_matrix;
}

//Fun�ao auxiliar: verifica se tr�s vetores formam uma base ortonormal
bool isOrthonormal(const glm::vec4& u, const glm::vec4& v, const glm::vec4& n) {
	if (dotproduct(u, v) != 0.0f || dotproduct(u, n) != 0.0f || dotproduct(v, n) != 0.0f) {
		return false;
	}
	if (norm(u) != 1.0f || norm(v) != 1.0f || norm(n) != 1.0f) {
		return false;
	}
	return true;
}
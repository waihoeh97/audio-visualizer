#ifndef _AUDIO_VISUALIZER_H
#define _AUDIO_VISUALIZER_H
#define PI 3.14159263
#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))
#define DEF_D 5
#define RECT_VERTICE_W 50
#define RECT_VERTICE_H 50

#include "demo_base.h"
#include "bitmap.h"
#include <cmath>
#include <vector>
#include <fmod.hpp>
#include <fmod_errors.h>

#define TEXTURE_COUNT 1

const int RECT_VERTEX_COUNT = (RECT_VERTICE_W - 1) * (RECT_VERTICE_H - 1) * 6;
const int RECT_VERTEX_ARRAY_SIZE = RECT_VERTEX_COUNT * 3;

#define SPECTRUM_SIZE 512

void ERRCHECK(FMOD_RESULT result)
{
	if (result != FMOD_OK)
	{
		//printf("FMOD error! (%d) %s\n", result, FMOD_ErrorString(result));
		std::cout << "FMOD ERROR" << FMOD_ErrorString(result) << std::endl;
	}
}

class Vertex
{
public:
	float x;
	float y;
	float z;

	Vertex()
	{
		x = 0;
		y = 0;
		z = 0;
	}

	Vertex(float a, float b, float c)
	{
		x = a;
		y = b;
		z = c;
	}
};

Vertex SphereFunction(float u, float v, float radius)
{
	u = u * PI / 180.0f;
	v = v * PI / 180.0f;
	Vertex answer;
	answer.x = cos(u) * sin(v) * radius;
	answer.y = cos(v) * radius;
	answer.z = sin(u) * sin(v) * radius;
	return answer;
}

class AudioVisualizer : public DemoBase, Vertex
{
private:

	FMOD::System* m_fmodSystem;
	FMOD::Sound* m_music;
	FMOD::Channel* m_musicChannel;
	float m_spectrumLeft[SPECTRUM_SIZE];
	float m_spectrumRight[SPECTRUM_SIZE];


	float m_rotation;
	float m_value;
	float m_translate;
	float m_frequency;

	float spectrum;

	float offSet;

	GLfloat mRectVertices[RECT_VERTEX_ARRAY_SIZE];

	GLuint mTextureID[TEXTURE_COUNT];

	void loadTexture(const char* path, GLuint textureID)
	{
		CBitmap bitmap(path);

		// Create Linear Filtered Texture
		glBindTexture(GL_TEXTURE_2D, textureID);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);	// apply texture wrapping
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// apply texture wrapping

																		//// nearest (minecraft) filtering
																		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); // near filtering
																		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // far filtering


																		// bilinear filtering
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // near filtering
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // far filtering


		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, bitmap.GetWidth(), bitmap.GetHeight(),
			0, GL_RGBA, GL_UNSIGNED_BYTE, bitmap.GetBits());
	}
public:
	float x;
	float y;
	float z;

	void init()
	{
		m_rotation += 0.05f;

		//=============init fmod==========================================================================
		FMOD_RESULT  result;
		unsigned int version;

		result = FMOD::System_Create(&m_fmodSystem);
		ERRCHECK(result);

		/*result = m_fmodSystem->getVersion(&version);
		ERRCHECK(result);

		if (version < FMOD_VERSION) {
		printf("FMOD Error! You are using and old version of FMOD.", version, FMOD_VERSION);
		}*/

		// Initialize fmod system
		result = m_fmodSystem->init(32, FMOD_INIT_NORMAL, 0);
		ERRCHECK(result);

		// load and set up music
		result = m_fmodSystem->createStream("../media/FuwaFuwaTime.mp3", FMOD_SOFTWARE, 0, &m_music);
		ERRCHECK(result);

		// play loaded mp3 music
		result = m_fmodSystem->playSound(FMOD_CHANNEL_FREE, m_music, false, &m_musicChannel);
		ERRCHECK(result);
		//=================================================================================================

		glGenTextures(TEXTURE_COUNT, mTextureID);
		loadTexture("../media/waterTexture.bmp", mTextureID[0]);
	}

	void deinit()
	{
	}

	void drawCube(float m_frequency)
	{
		glBegin(GL_TRIANGLES);

		//glColor3f(255.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, m_frequency, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, m_frequency, 0.0f);
		glVertex3f(1.0f, m_frequency, 0.0f);

		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(1.0f, m_frequency, 1.0f);
		glVertex3f(1.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, m_frequency, 1.0f);
		glVertex3f(1.0f, m_frequency, 1.0f);

		// x-axis
		//glColor3f(255.0f, 69.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, m_frequency, 1.0f);
		glVertex3f(0.0f, m_frequency, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, m_frequency, 1.0f);


		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, m_frequency, 1.0f);
		glVertex3f(1.0f, m_frequency, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 1.0f);
		glVertex3f(1.0f, m_frequency, 1.0f);

		// y-axis
		//glColor3f(255.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, m_frequency, 0.0f);
		glVertex3f(1.0f, m_frequency, 0.0f);
		glVertex3f(1.0f, m_frequency, 1.0f);
		glVertex3f(0.0f, m_frequency, 0.0f);
		glVertex3f(0.0f, m_frequency, 1.0f);
		glVertex3f(1.0f, m_frequency, 1.0f);

		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 1.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);

		glEnd();
	}

	void drawCircle(float k, float r, float h)
	{
		glBegin(GL_LINES);
		glColor3f(128.0f, 0.0f, 128.0f);
		for (int i = 0; i < 180; i++)
		{
			x = r * cos(i) - h;
			y = r * sin(i) + k;
			glVertex3f(x + k, y - h, 0);

			x = r * cos(i + 0.1) - h;
			y = r * sin(i + 0.1) + k;
			glVertex3f(x + k, y - h, 0);
		}
		glEnd();
	}

	void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius) {
		int i;
		int triangleAmount = 20;

		GLfloat twicePi = 2.0f * M_PI;

		glBegin(GL_TRIANGLE_FAN);
		glColor3f(128.0f, 0.0f, 128.0f);
		glVertex2f(x, y); // center of circle
		for (i = 0; i <= triangleAmount; i++) {
			glVertex2f(
				x + (radius * cos(i *  twicePi / triangleAmount)),
				y + (radius * sin(i * twicePi / triangleAmount))
			);
		}
		glEnd();
	}

	void DrawSphere(float startU, float startV, float endU, float endV, float radius, float UResolution, float VResolution, float amplifier)
	{
		float stepU = (endU - startU) / UResolution;
		float stepV = (endV - startV) / VResolution;

		glBegin(GL_TRIANGLES);
		for (int i = 0; i < UResolution; i++)
		{
			for (int j = 0; j < VResolution; j++)
			{
				float u = i * stepU + startU;
				float v = j * stepV + startV;
				float un = (i + 1 == UResolution) ? endU : (i + 1)*stepU + startU;
				float vn = (j + 1 == VResolution) ? endV : (j + 1)*stepV + startV;

				Vertex p0 = SphereFunction(u, v, radius + amplifier);
				Vertex p1 = SphereFunction(u, vn, radius + amplifier);
				Vertex p2 = SphereFunction(un, v, radius + amplifier);
				Vertex p3 = SphereFunction(un, vn, radius + amplifier);

				glVertex3f(p0.x, p0.y, p0.z);
				glVertex3f(p2.x, p2.y, p2.z);
				glVertex3f(p1.x, p1.y, p1.z);

				glVertex3f(p3.x, p3.y, p3.z);
				glVertex3f(p1.x, p1.y, p1.z);
				glVertex3f(p2.x, p2.y, p2.z);
			}
		}
		glEnd();
	}

	void genWave(float offsetRadian)
	{
		const float width = 5.0f;
		const float height = 5.0f;
		const float halfWidth = width * 0.5f;
		const float halfHeight = height * 0.5f;
		const float waveAmpMul = 0.1f;
		const float xToRadianMul = 5.0f;

		m_translate += 0.01f;

		int currentVert = -1;
		int currentIndex = -1;

		offsetRadian += 1.0f;

		for (int h = 0; h<RECT_VERTICE_H - 1; h++)
		{
			for (int w = 0; w<RECT_VERTICE_W - 1; w++)
			{
				float entityRotation = m_translate - (0.1f * (float)h);

				//========= 6 vertices to form one sub-rectangle

				//1st vertex
				int vertex1 = ++currentVert;
				mRectVertices[vertex1 * 3] = -halfWidth + (float)(w) / (float)RECT_VERTICE_W * width;
				mRectVertices[vertex1 * 3 + 1] = waveAmpMul * cos(mRectVertices[vertex1 * 3] * xToRadianMul + offsetRadian);
				mRectVertices[vertex1 * 3 + 2] = -halfHeight + (float)(h) / (float)RECT_VERTICE_H * height;

				//2nd vertex
				int vertex2 = ++currentVert;
				mRectVertices[vertex2 * 3] = -halfWidth + (float)(w) / (float)RECT_VERTICE_W * width;
				mRectVertices[vertex2 * 3 + 1] = waveAmpMul * cos(mRectVertices[vertex2 * 3] * xToRadianMul + offsetRadian);
				mRectVertices[vertex2 * 3 + 2] = -halfHeight + (float)(h + 1) / (float)RECT_VERTICE_H * height;

				//3rd vertex
				int vertex3 = ++currentVert;
				mRectVertices[vertex3 * 3] = -halfWidth + (float)(w + 1) / (float)RECT_VERTICE_W * width;
				mRectVertices[vertex3 * 3 + 1] = waveAmpMul * cos(mRectVertices[vertex3 * 3] * xToRadianMul + offsetRadian);
				mRectVertices[vertex3 * 3 + 2] = -halfHeight + (float)(h + 1) / (float)RECT_VERTICE_H * height;

				//4th vertex
				int vertex4 = ++currentVert;
				mRectVertices[vertex4 * 3] = mRectVertices[vertex3 * 3];
				mRectVertices[vertex4 * 3 + 1] = mRectVertices[vertex3 * 3 + 1];
				mRectVertices[vertex4 * 3 + 2] = mRectVertices[vertex3 * 3 + 2];

				//5th vertex
				int vertex5 = ++currentVert;
				mRectVertices[vertex5 * 3] = -halfWidth + (float)(w + 1) / (float)RECT_VERTICE_W * width;
				mRectVertices[vertex5 * 3 + 1] = waveAmpMul * cos(mRectVertices[vertex5 * 3] * xToRadianMul + offsetRadian);
				mRectVertices[vertex5 * 3 + 2] = -halfHeight + (float)(h) / (float)RECT_VERTICE_H * height;

				//6th vertex
				int vertex6 = ++currentVert;
				mRectVertices[vertex6 * 3] = mRectVertices[vertex1 * 3];
				mRectVertices[vertex6 * 3 + 1] = mRectVertices[vertex1 * 3 + 1];
				mRectVertices[vertex6 * 3 + 2] = mRectVertices[vertex1 * 3 + 2];
			}
		}
	}

	void drawWorldOriginLines(const Matrix& viewMatrix)
	{
		// draw axis.
		glLoadMatrixf((GLfloat*)viewMatrix.mVal);
		glBegin(GL_LINES);
		glColor3f(1.0f, 0.3f, 0.3f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);

		glColor3f(0.3f, 1.0f, 0.3f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);

		glColor3f(0.3f, 0.3f, 1.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glEnd();
	}

	void drawCubeHuge(void)
	{
		offSet += 0.01f;
		// Prevent buffer overflow
		if (offSet > 2.0f * PI)
		{
			offSet -= 2.0f * PI;
		}
		genWave(offSet);

		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		//glEnable(GL_BLEND);
		//glEnable(GL_TEXTURE_2D);
		//glEnable(GL_LIGHTING);

		glBindTexture(GL_TEXTURE_2D, mTextureID[3]);

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, mRectVertices);

		glDrawArrays(GL_TRIANGLES, 0, RECT_VERTEX_COUNT);

		glDisableClientState(GL_VERTEX_ARRAY);
	}

	void updateFmod(float m_frequency, const Matrix& viewMatrix)
	{
		m_fmodSystem->update();
		// Get spectrum for left and right stereao channels
		m_musicChannel->getSpectrum(m_spectrumLeft, SPECTRUM_SIZE, 0, FMOD_DSP_FFT_WINDOW_RECT);
		m_musicChannel->getSpectrum(m_spectrumLeft, SPECTRUM_SIZE, 1, FMOD_DSP_FFT_WINDOW_RECT);

		//this is the frequency
		// Right
		for (int i = 0; i < 5; i++)
		{
			m_frequency = m_spectrumLeft[i] + m_spectrumRight[i] * 0.5f;
			if (m_frequency >= 0.25f)
			{
				glColor3f(255.0f, 255.0f, 0.0f);
			}
			else
			{
				glColor3f(255.0f, 0.0f, 0.0f);
			}
			Matrix translate = Matrix::makeTranslationMatrix(i + 1.0f, 0.0f, 0.5f);
			Matrix scale = Matrix::makeScaleMatrix(1.0f, 3.0f, 1.0f);
			Matrix viewSpaceMatrix = viewMatrix * translate * scale;
			glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
			drawCube(m_frequency);
		}

		// Left
		for (int i = 0; i < 5; i++)
		{
			m_frequency = m_spectrumLeft[i] + m_spectrumRight[i] * 0.5f;
			if (m_frequency >= 0.25f)
			{
				glColor3f(255.0f, 255.0f, 0.0f);
			}
			else
			{
				glColor3f(255.0f, 0.0f, 0.0f);
			}
			Matrix translate = Matrix::makeTranslationMatrix(i - 6.0f, 0.0f, 0.5f);
			Matrix scale = Matrix::makeScaleMatrix(1.0f, 3.0f, 1.0f);
			Matrix viewSpaceMatrix = viewMatrix * translate * scale;
			glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
			drawCube(m_frequency);
		}

		// Middle Sphere
		m_frequency = m_spectrumLeft[1] + m_spectrumRight[1] * 0.5f;
		if (m_frequency >= 0.25f)
		{
			glColor3f(230.0f, 230.0f, 250.0f);
		}
		else
		{
			glColor3f(128.0f, 0.0f, 128.0f);
		}
		Matrix sphere = Matrix::makeScaleMatrix(1.0f, 1.0f, 1.0f) * Matrix::makeTranslationMatrix(0.0f, 0.0f, 1.0f);
		Matrix viewSpaceMatrix = viewMatrix * sphere;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		//drawFilledCircle(0.0f, 0.0f, m_frequency);
		DrawSphere(0.0f, 0.0f, 360.0f, 180.0f, 1.0f, 20.0f, 20.0f, m_frequency);

		// Wave bottom
		glColor3f(0.0f, 191.0f, 255.0f);
		offSet += m_frequency;
		Matrix wave = Matrix::makeScaleMatrix(5.0f, 5.0f, 6.0f) * Matrix::makeTranslationMatrix(0.0f, -0.5f, 0.0f);
		viewSpaceMatrix = viewMatrix * wave;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		drawCubeHuge();

		// Wave top
		glColor3f(0.0f, 191.0f, 255.0f);
		offSet += m_frequency;
		wave = Matrix::makeScaleMatrix(5.0f, 5.0f, 6.0f) * Matrix::makeTranslationMatrix(0.0f, 1.5f, 0.0f);
		viewSpaceMatrix = viewMatrix * wave;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		drawCubeHuge();

		// Wave left
		glColor3f(0.0f, 0.0f, 255.0f);
		wave = Matrix::makeScaleMatrix(5.0f, 5.0f, 6.5f) * Matrix::makeRotateMatrix(90.0f, Vector(0.0f, 0.0f, 1.0f)) * Matrix::makeTranslationMatrix(0.0f, -2.0f, 0.0f);
		viewSpaceMatrix = viewMatrix * wave;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		drawCubeHuge();

		// Wave right
		glColor3f(0.0f, 0.0f, 255.0f);
		wave = Matrix::makeScaleMatrix(5.0f, 5.0f, 6.5f) * Matrix::makeRotateMatrix(90.0f, Vector(0.0f, 0.0f, 1.0f)) * Matrix::makeTranslationMatrix(0.0f, 2.0f, 0.0f);
		viewSpaceMatrix = viewMatrix * wave;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		drawCubeHuge();

		// Wave front
		glColor3f(0.0f, 191.0f, 255.0f);
		offSet += m_frequency;
		wave = Matrix::makeScaleMatrix(5.0f, 5.0f, 5.0f) * Matrix::makeRotateMatrix(90.0f, Vector(1.0f, 0.0f, 0.0f)) * Matrix::makeTranslationMatrix(0.0f, 2.0f, 0.0f);
		viewSpaceMatrix = viewMatrix * wave;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		drawCubeHuge();

		// Wave back
		glColor3f(0.0f, 191.0f, 255.0f);
		offSet += m_frequency;
		wave = Matrix::makeScaleMatrix(5.0f, 5.0f, 5.0f) * Matrix::makeRotateMatrix(90.0f, Vector(1.0f, 0.0f, 0.0f)) * Matrix::makeTranslationMatrix(0.0f, -3.0f, 0.0f);
		viewSpaceMatrix = viewMatrix * wave;
		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		drawCubeHuge();
	}

	void draw(const Matrix& viewMatrix)
	{
		updateFmod(0, viewMatrix);

		// draw axis.
		drawWorldOriginLines(viewMatrix);

		Matrix rotate1 = Matrix::makeRotateMatrix(0, Vector(0.0f, 0.0f, 1.0f));
		Matrix translate1 = Matrix::makeTranslationMatrix(0.0f, 0.0f, 0.0f);

		// NOTE on OpenGL Matrix model:
		// Screen = Proj * View * Model
		// Model = TransformA(3rd) * TransformB(2nd) * TransformC(1st) (Trasnform could be Rotate, Scale, Translate, etc)

		// perform model transformation
		Matrix modelMatrix = translate1 * rotate1;

		// convert model space to view space
		Matrix viewSpaceMatrix = viewMatrix * modelMatrix;

		glLoadMatrixf((GLfloat*)viewSpaceMatrix.mVal);
		//drawCubeHuge();
	}
};

#endif

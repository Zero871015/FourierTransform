#include "DataManager.h"

DataManager::DataManager(int h, int w)
{
	ImageHeight = h;
	ImageWidth = w;

	InputImage = new int*[h];
	OutputImage = new int*[h];
	FreqReal = new double*[h];
	FreqImag = new double*[h];

	for (int i = 0; i < h; i++)
	{
		InputImage[i] = new int[w];
		OutputImage[i] = new int[w];
		FreqReal[i] = new double[w];
		FreqImag[i] = new double[w];

	}

	for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++)
		{
			InputImage[i][j] = 0;
			OutputImage[i][j] = 0;
			FreqReal[i][j] = 0;
			FreqImag[i][j] = 0;
		}
}

void DataManager::SetPixel(int x, int y, int pixelValue)
{
	InputImage[y][x] = pixelValue;
}

void DataManager::SetFreqReal(int x, int y, double value)
{
	FreqReal[y][x] = value;
}

void DataManager::SetFreqImag(int x, int y, double value)
{
	FreqReal[y][x] = value;
}

int DataManager::GetImageHeight()
{
	return ImageHeight;
}

int DataManager::GetImageWidth()
{
	return ImageWidth;
}

int ** DataManager::GetInputImage()
{
	return InputImage;
}

int ** DataManager::GetOutputImage()
{
	return OutputImage;
}

double ** DataManager::GetFreqReal()
{
	return FreqReal;
}

double ** DataManager::GetFreqImag()
{
	return FreqImag;
}

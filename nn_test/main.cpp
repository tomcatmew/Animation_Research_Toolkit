#include <GLFW/glfw3.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <stdlib.h>
#include <deque>
#include <torch/torch.h>
#include <torch/script.h>

// my utilities include folder files
#include "ycmatrix.cpp"
#include "ycdraw.cpp"

using namespace YCutility;
using namespace YCdraw;


int main(void)
{
	std::cout << std::fixed;
	torch::jit::script::Module module;
	try {
		module = torch::jit::load("traced.pt");
	}
	catch (const c10::Error& e) {
		std::cerr << "error load model \n";
		return -1;
	}

	std::vector<torch::jit::IValue> inputs;

	// ==== read the ininput 
	std::fstream in("test.bvh");
	std::string line;
	std::vector<float> float_vector;

	std::getline(in, line);
	float readvalue;
	std::stringstream ss(line);

	while (ss >> readvalue)
	{
		float_vector.push_back(readvalue);
	}

	std::cout << std::endl << "input size:" << float_vector.size() << std::endl;
	std::cout << "input :" << std::endl;
	for (std::vector<float>::const_iterator i = float_vector.begin(); i != float_vector.end(); ++i)
		std::cout << *i << ' ';

	// ====

	torch::Tensor input_para = torch::from_blob(float_vector.data(), { 96 });

	//std::cout << "print Input tensor" << std::endl;
	//std::cout << std::fixed << input_para << std::endl;
	inputs.push_back(input_para);

	at::Tensor output_para = module.forward(inputs).toTensor();
	std::cout << "print Output tensor size: " << std::endl;

	//std::cout << output_para.sizes()[0] << std::endl;

	std::vector<double> last_vector;
	for (int i = 0; i < output_para.sizes()[0]; i++)
	{
		double pos_val = output_para[i].item<double>();

		//if (i == 1) {
		//	last_vector.push_back(pos_val - 180.0f);
		//}
		//else {
			last_vector.push_back(pos_val);
		//}
	}
	for (std::vector<double>::const_iterator i = last_vector.begin(); i != last_vector.end(); ++i)
		std::cout << *i << ' ';


	for (int j = 0; j < 20; j++)
	{
		float_vector.clear();
		for (int i = 0; i < last_vector.size(); i++) {
			float_vector.push_back((float)last_vector[i]);
		}
		torch::Tensor input_tens = torch::from_blob(float_vector.data(), { 96 });

		inputs.clear();
		inputs.push_back(input_tens);

		at::Tensor output_tens = module.forward(inputs).toTensor();

		last_vector.clear();

		for (int i = 0; i < output_tens.sizes()[0]; i++)
		{
			double new_pos_val = output_tens[i].item<double>();

			//if (i == 1) {
			//	last_vector.push_back(new_pos_val - 180.0f);
			//}
			//else {
				last_vector.push_back(new_pos_val);
			//}
		}

		int tempt_count = 0;
		for (std::vector<double>::const_iterator i = last_vector.begin(); i != last_vector.end(); ++i)
		{
			std::cout << *i << ' ';
			tempt_count++;
			if (tempt_count > 10) 
				break;
		}
		std::cout << std::endl;
	}
}
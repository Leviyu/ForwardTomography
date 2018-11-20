#include "../02_forward_tomo/forward_tomography.h"
#include "../01_cpp_lib/hongyulibcpp.h"

using namespace std;
int main()
{
	string PHASE;
	cout << "============= forward tomography begin ========================== " << endl;

	string record_file;
	string model_path;
	string input_file;
	string MODEL;
	
	ifstream myfile;
	string infile = "input";
	myfile.open(infile.c_str());

	myfile >> record_file >> model_path >> input_file >> MODEL;

	new_tomo my_tomo;
	my_tomo.INFILE = input_file;
	my_tomo.read_INFILE();
	my_tomo.flag_is_new_tomo = "new";
	my_tomo.read_standard_tomo();


	// define increment 
	my_tomo.deg_inc = 10;
	//double deg_inc = 20;
	//double deg_inc = 30;


	big_new_record my_big_record;
	my_tomo.current_iteration = 4;
	my_tomo.get_sig_index();
	my_tomo.store_smoother(&my_big_record);
	my_tomo.smooth_dvs_model(&my_big_record);

	my_tomo.output_starting_tomography();




	//my_tomo.construct_RMS_profile("on");
	//my_tomo.output_RMS_profile();


	// Read in record
	//my_big_record.sta_num = count_file_num(record_file);
	//my_big_record.initiate_big_record();
	//my_big_record.my_record.resize(my_big_record.sta_num);
	//my_big_record.record_file = record_file;
	//my_big_record.read_record_file(my_big_record.my_record);
	////my_big_record.big_record_read_cross_point_file(&my_tomo);

	// output starting model
	//my_tomo.output_tomography_info3("start");
	//my_tomo.forward_tomography_func(&my_big_record);
	//my_tomo.output_tomography_info3("final.3");



	// OUTPUT
	//my_tomo2.output_starting_tomography();
	//my_tomo.output_tomography_info2();
	//my_tomo.output_time_info();


	//my_tomo.output_record_path_crosssection(&my_big_record);
	//my_tomo.output_record_path_line(&my_big_record);
	//my_tomo.output_record_path_line_in_one_file(&my_big_record);

	cout << "============= forward tomography ends  ========================== " << endl;
	return 0;
}

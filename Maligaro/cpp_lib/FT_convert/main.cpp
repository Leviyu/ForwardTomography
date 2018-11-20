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
	

	cin >> record_file >> model_path >> input_file >> MODEL;

	new_tomo tomo;
	new_tomo my_tomo;
	new_tomo my_tomo2;

	tomo.INFILE = input_file;
	tomo.read_INFILE();
	tomo.flag_is_new_tomo = "old";
	tomo.MODEL_INFILE = model_path + "/" + MODEL + ".INFILE";
	tomo.read_MODEL_INFILE();
	tomo.initiate_tomo();
	tomo.MODEL_model = model_path + "/" + MODEL + ".model";
	tomo.read_tomography();


	// convert tomo into standard tomography
	my_tomo.flag_is_new_tomo = "new";
	//my_tomo2.flag_is_new_tomo = "new";

	my_tomo.convert_to_new_tomo(&tomo);
	//my_tomo2.convert_to_new_tomo(&tomo);

	//my_tomo.construct_RMS_profile();

	// Read in record
	//big_new_record my_big_record;
	//my_big_record.sta_num = count_file_num(record_file);

	//my_big_record.initiate_big_record();
	//my_big_record.record_file = record_file();
	//my_big_record.read_record_file();
	//my_big_record.big_record_read_cross_point_file(&my_tomo);

	//my_tomo.forward_tomography_func(&my_big_record);



	// OUTPUT
	//my_tomo.output_delta_tomography(&my_tomo2);
	//my_tomo2.output_starting_tomography();
	//my_tomo.output_starting_tomography();
	my_tomo.output_tomography_info_standard();
	//my_tomo.output_time_info();


	//my_tomo.output_record_path_crosssection(&my_big_record);
	//my_tomo.output_record_path_line(&my_big_record);
	//my_tomo.output_record_path_line_in_one_file(&my_big_record);

	cout << "============= forward tomography ends  ========================== " << endl;
	return 0;
}

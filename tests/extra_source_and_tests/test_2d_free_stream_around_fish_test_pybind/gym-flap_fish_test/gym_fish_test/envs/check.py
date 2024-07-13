# import test_2d_flow_stream_around_fish_test_pybind as test_fish
import test_2d_flow_stream_around_fish_flap_pybind as train
import test_2d_flow_stream_around_fish_flap_test_pybind as test

b = test.from_sph_relaxation_and_test(101)  # generate configuration xml fish file
b = test.from_sph_reload_and_test(101)



for i in range(32):
	pressure = b.GetFishPressurePoint(i)
	v_x = b.GetFishVelocityX(i)
	v_y = b.GetFishVelocityY(i)
	pos_x = b.GetFishPositionX(i)
	pos_y = b.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)
b.RunCase(101, 0.5)
for i in range(32):
	pressure = b.GetFishPressurePoint(i)
	v_x = b.GetFishVelocityX(i)
	v_y = b.GetFishVelocityY(i)
	pos_x = b.GetFishPositionX(i)
	pos_y = b.GetFishPositionY(i)
	print(i, ' pressurePoint: ', pressure, ' V_x: ', v_x, ' V_y: ', v_y, ' pos_x: ', pos_x, ' pos_y: ', pos_y)






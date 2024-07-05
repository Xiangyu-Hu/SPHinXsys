# import test_2d_flow_stream_around_fish_test_pybind as test_fish
import test_2d_flow_stream_around_fish_flap_pybind as train
import test_2d_flow_stream_around_fish_flap_test_pybind as test

r = test.from_sph_relaxation_and_test(1)
a = test.from_sph_reload_and_test(1)

# r = train.from_sph_relaxation(1)
# a = train.from_sph_reload_and_train(1)

# b = test.from_sph_reload_and_test(101)
#
# a.SetFreq(1) # 2.5 - 2.7 between 5 times get average time of each parameter
# 2.5 outside fish speed < fluid speed last 6.1S
# 2.8  fish speed < fluid speed last 1.3S each situation different, need to change to different 2nd:
# 3  fish speed > fluid speed
# a.SetLambda(1)
a.SetAmUp(0.12)
a.SetAmDown(0.12)

a.RunCase(1, 20)

# action_time = 0
# #
# for i in range(16):
# 	pressure_ = a.GetFishPressurePoint(i)
# 	p_x = a.GetFishPositionX(i)
# 	p_y = a.GetFishPositionY(i)
# 	v_x = a.GetFishVelocityX(i)
# 	v_y = a.GetFishVelocityY(i)

# for times in range(1500):
#     for i in range(30):
#         action_time += 0.05 * 0.1 / 3
#         a.SetAmUp(0.12)
#         a.SetAmDown(0.12)
#         a.RunCase(1, action_time)
#         # for pt in range(16):
#         #     pressure_ = a.GetFishPressurePoint(pt)
#         #     p_x = a.GetFishPositionX(pt)
#         #     p_y = a.GetFishPositionY(pt)
#         #     v_x = a.GetFishVelocityX(pt)
#         #     v_y = a.GetFishVelocityY(pt)
#         #     strain = a.GetFishActiveStrain(pt)
#         #     print(i, ' pressure_: ', pressure_, ' p_x: ', p_x, ' p_y: ', p_y, ' v_x: ', v_x, ' v_y: ', v_y, ' strain: ', strain)
#         v_x = a.GetFishVelocityX(7)
#         v_y = a.GetFishVelocityY(7)
#         pressure_ = a.GetFishPressurePoint(7)
#         p_x = a.GetFishPositionX(7)
#         p_y = a.GetFishPositionY(7)
#         strain = a.GetFishActiveStrain(7)
#         # print( times* 0.05 + 0.05 * 0.1* i , ' V_x: ', v_x, '   V_y: ', v_y)
#         file = open(f'fish_v.txt', 'a')
#         file.write('  ')
#         file.write(str(round(times* 0.05 + 0.05 * 0.1* i, 4)))
#         file.write('  , ')
#         file.write(str(v_x))
#         file.write('  ,  ')
#         file.write(str(v_y))
#         file.write('  ,  ')
#         file.write(str(p_x))
#         file.write('  ,  ')
#         file.write(str(p_y))
#         file.write('  ,  ')
#         file.write(str(pressure_))
#         file.write('  ,  ')
#         file.write(str(strain))
#         file.write('\n')

# b.RunCase(101, 0.3)
#  1 1.5 0.15 0.15   3 left
#  1  1  0.15 0.15  1 left 2 down  2 upside
#  1  1  0.12 0.12  3 left 1 down  2 upside  fast
#  1.5 0.5 0.12 0.12  3 down 1 left 1 right center 1
#  1.5 0.5 0.15 0.15  4 left 1 down
#  1.5 0.5  0.1 0.1   1 right  2 left 1 down


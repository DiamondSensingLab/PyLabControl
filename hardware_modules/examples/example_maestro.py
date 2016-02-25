import time
import sys
sys.path.insert(0,"C:\\Users\Experiment\\PycharmProjects\\PythonLab\\")
import hardware_modules.maestro as maestro

# define com port to communicate with servo (check in device manager for pololu  controler command port)
servo = maestro.Controller('COM5')

# define what to test

test_case = 'filterwheel' # beamblock, motor, controller, filterwheel


# set channel
channel = 0

# ================== test filter wheel =======================

if test_case == 'filterwheel':
    filter = maestro.FilterWheel(servo, channel, {'ND1.0': 4*600, 'LP':4*1550, 'ND2.0':4*2500})
    filter.goto('ND2.0')
    # filter.goto('ND1.0')
    #filter.goto('LP')
    # block1.block()
    # close communication channel
    # servo.close()

    for i in range(10):
        print(servo.getPosition(channel))

# ================== test controler =======================
# # set ranges, acceleration etc.
if test_case == 'controller':
    # servo.setRange(channel, 2000, 12000)
    # servo.setAccel(channel, 0)
    # servo.setSpeed(channel, 0)
    # servo.setTarget(channel, 6100)
    # servo.setTarget(channel, 5000)
    servo.disable(channel)

# ================== test motor =======================
if test_case == 'motor':
    motor = maestro.Motor(servo, channel)
    motor.rotate(1000)

    time.sleep(2)
    # motor.rotate(0)

    motor.stop()
# ================== test beam block =======================
if test_case == 'beamblock':
    channel = 4
    block1 = maestro.BeamBlock(servo, channel)
    # block1.open()
    block1.block()

    # for i in range(10):
    #     print(servo.getPosition(channel))
    # time.sleep(0.2)
    servo.setTarget(0,  channel)
    # servo.goHome()
    # close communication channel
    # servo.close()



import tensorflow as tf

print('gpu detected', tf.test.is_gpu_available(), flush=True)
print('gpu current device', tf.config.list_physical_devices('GPU'), flush=True)
import torch

print('gpu detected', torch.cuda.is_available(), flush=True)
print('gpu device count', torch.cuda.device_count(), flush=True)
print('gpu current device', torch.cuda.current_device(), flush=True)
print('gpu device name', torch.cuda.get_device_name(0), flush=True)
from diffusers import DiffusionPipeline, StableDiffusionXLImg2ImgPipeline
import torch

prj_path = "username/repo_name"
model = "stabilityai/stable-diffusion-xl-base-1.0"
pipe = DiffusionPipeline.from_pretrained(model,torch_dtype=torch.float16,)
pipe.to("cuda")
pipe.load_lora_weights("pytorch_lora_weights.safetensors")

refiner = StableDiffusionXLImg2ImgPipeline.from_pretrained("stabilityai/stable-diffusion-xl-refiner-1.0",torch_dtype=torch.float16)
refiner.to("cuda")

prompt = "photo of sameer in suit for CV"

seed = 42
generator = torch.Generator("cuda").manual_seed(seed)
image = pipe(prompt=prompt, generator=generator).images[0]
image = refiner(prompt=prompt, generator=generator, image=image).images[0]
image.save(f"generated_image.png")
from cv2 import cv2
import os

image_folder = 'C:\\Users\\felix\\Documents\\GitHub\\gymnasiearbete\\images\\070520-041042'
video_name = 'C:\\Users\\felix\\Documents\\GitHub\\gymnasiearbete\\videos\\070520-041042.avi'

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

video = cv2.VideoWriter(video_name, 0, 5, (width,height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
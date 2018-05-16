from PIL import Image, ImageDraw, ImageFont

background = Image.new('RGB', color = (255, 255, 255), size = (2000, 1000))
add_text = ImageDraw.Draw(background)

# Andale Mono
font = ImageFont.truetype("/Library/Fonts/Arial Bold.ttf", 18)
hgts = ['L12-->B2', 'C1-->n22', 'C1-->D5', 'L6-->n11', 'L6-->G24', 'G8-->G7', 'G5-->G8', 'G5-->G26', 'G21-->G23', 'n19-->G21', 'G16-->G25', 'L12-->B2', 'C1-->n22', 'C1-->D5', 'L6-->n11', 'L6-->G24', 'G8-->G7', 'G5-->G8', 'G5-->G26', 'G21-->G23', 'n19-->G21', 'G16-->G25']
x_start = 30
y_start = 80

for hgt in hgts:

    if x_start < 900:
        if hgt == 'C1-->D5':
            add_text.text((x_start, y_start), hgt, (225, 0, 0), font = font)
            x_start += 120
        else:
            add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font)
            x_start += 120
    elif x_start > 900:
        x_start = 30
        y_start += 50
        if hgt == 'C1-->D5':
            add_text.text((x_start, y_start), hgt, (225, 0, 0), font = font)
            x_start += 120
        else:
            add_text.text((x_start, y_start), hgt, (0, 0, 0), font = font)
            x_start += 120

background.save("/Users/weizhisong/Desktop/a_test.png")

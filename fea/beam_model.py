from cadgen import build123d as bd
from cadgen import step

L, b, h = 200.0, 30.0, 6.0


@step(out="beam.step")
def beam():
    body = bd.Box(L, b, h)
    body.label = "beam"
    return body


if __name__ == "__main__":
    beam()

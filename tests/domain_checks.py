from datapack.domain import transform,compose,signed_area,segment_length,interpolate
def check(case):
    before=[.001,0,0,0,.001,0];after=[-1,0,20,0,1,0]
    case.assertEqual(transform(compose(after,before),[1000,2000]),[19,2])
    case.assertEqual(signed_area([[0,0],[1,0],[1,1],[0,1]]),1)
    case.assertEqual(segment_length([0,0],[3,4]),5)
    case.assertEqual(interpolate([0,10],[10,20],.5),[5,15])
    with case.assertRaises(ValueError):interpolate([0],[1],2)

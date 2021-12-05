#!/usr/bin/env python


from argparse import ArgumentParser

from scipy.integrate import solve_ivp


class Model:
    """Generic Model base class"""
    def __init__(self, **kwargs):
        params = dict(self.parameters)
        ics = dict(self.variables)
        for arg, val in kwargs.items():
            if arg in params:
                params[arg] = value
            elif arg in ics:
                ics[arg] = value
            else:
                raise ValueError(f"No such argument {arg}")

        self.X0 = tuple(ics[name] for name, _ in self.variables)
        self.P = tuple(params[name] for name, _ in self.parameters)

        for name, _ in self.parameters:
            setattr(self, name, params[name])
        for index, (name, _) in enumerate(self.variables):
            setattr(self, name, index)

    @classmethod
    def main(cls, args):
        parser = ArgumentParser()
        parser.add_argument('--filename')
        for name, default in cls.parameters + cls.variables:
            parser.add_argument(f'--{name}', default=default, type=float)

        args = parser.parse_args(args)
        args = vars(args) # to dict

        filename, kwargs = args.pop('filename'), args
        model = cls(**kwargs)

        t = np.linspace(0, 10, 100)

    def f_tX(self, t, X):
        # format expected for solve_ivp
        dXdt = zeros_like(X)
        self.f_rhs(X, t, 

    def f_rhs(self, t, X, dXdt):
        raise NotImplementedError


class SHM(Model):
    """Simple harmonic motion model.

    x'' + b x' + w^2 x = 0

    State space form:

    x' = v
    v' = - b v - w^2 x
    """
    parameters = [
        ('b', 1),     # 1/s
        ('w', 1), # 1/s
    ]
    variables = [
        ('x', 1),     # m
        ('v', 0),     # m/s
    ]

    def f(self, X, t, dXdt):
        x, v = X
        b, w = self.P

        dxdt = v
        dvdt = - self.b*v - self.w**2*x

        dXdt[:] = [dxdt, dvdt]


if __name__ == "__main__":
    import sys
    SHM.main(sys.argv[1:])

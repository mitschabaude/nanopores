def check_solve(p):
    iterations = p.newton_solve()
    if iterations > 15:
        raise SystemExit("Error: Newton didn't converge in 15 iterations")



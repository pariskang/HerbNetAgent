"""
herbnetagent.cli — command-line entry point.

    python -m herbnetagent assess --codrug rivaroxaban --perpetrator-mg 0.006
    python -m herbnetagent assess --disease "venous thrombosis" --json out.json
"""
from __future__ import annotations
import argparse, json, sys
from .orchestrator import HerbNetAgent
from .knowledge import Knowledge


def _cmd_assess(args):
    agent = HerbNetAgent(Knowledge())
    res = agent.assess(disease=args.disease, codrug=args.codrug,
                       codrug_perpetrator_mg=args.perpetrator_mg,
                       use_network=not args.no_network)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(res, f, ensure_ascii=False, indent=1)
        print(f"wrote {args.json}")
    if args.md:
        with open(args.md, "w") as f:
            f.write(HerbNetAgent.render_markdown(res))
        print(f"wrote {args.md}")
    if not (args.json or args.md):
        print(HerbNetAgent.render_markdown(res))


def _cmd_info(args):
    agent = HerbNetAgent(Knowledge())
    print("HerbNetAgent 2.0 capabilities:", ", ".join(agent.reg.capabilities()))
    k = Knowledge()
    print(f"formula targets cached: {len(k.formula_targets())}")
    print(f"known co-drugs: {', '.join(k.KNOWN_CODRUGS)}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="herbnetagent",
                                description="Quantitative multi-scale TCM intelligence agent")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("assess", help="run an end-to-end assessment")
    a.add_argument("--disease", default="venous thrombosis")
    a.add_argument("--codrug", default=None, help="co-administered drug (e.g. rivaroxaban)")
    a.add_argument("--perpetrator-mg", type=float, default=None,
                   help="delivered perpetrator (e.g. quercetin) dose in mg for PBPK")
    a.add_argument("--no-network", action="store_true", help="skip network-proximity (offline)")
    a.add_argument("--json", default=None)
    a.add_argument("--md", default=None)
    a.set_defaults(func=_cmd_assess)

    i = sub.add_parser("info", help="show capabilities & cached knowledge")
    i.set_defaults(func=_cmd_info)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

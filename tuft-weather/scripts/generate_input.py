#!/usr/bin/env python3
import argparse
import random
from datetime import datetime, timedelta
from math import cos, pi

continent_locations = {
    "Europe": [
        ("Greece",        "Athens",        37.98),
        ("Germany",       "Berlin",        52.52),
        ("France",        "Paris",         48.86),
        ("Spain",         "Madrid",        40.42),
        ("Italy",         "Rome",          41.90),
        ("United Kingdom","London",        51.50),
    ],
    "Middle East": [
        ("Israel",        "Tel Aviv",      32.08),
        ("United Arab Emirates", "Dubai",  25.20),
        ("Saudi Arabia",  "Riyadh",        24.71),
        ("Iran",          "Tehran",        35.69),
    ],
    "Asia": [
        ("China",         "Beijing",       39.90),
        ("China",         "Chengdu",       30.67),
        ("China",         "Guangzhou",     23.13),
        ("China",         "Shanghai",      31.23),
        ("China",         "Shenyang",      41.80),
        ("Japan",         "Tokyo",         35.68),
        ("India",         "Mumbai",        19.08),
        ("South Korea",   "Seoul",         37.57),
        ("Thailand",      "Bangkok",       13.75),
        ("Malaysia",      "Kuala Lumpur",   3.14),
    ],
    "Africa": [
        ("Malawi",        "Lilongwe",     -13.96),
        ("Mozambique",    "Maputo",       -25.97),
        ("Namibia",       "Windhoek",     -22.56),
        ("Nigeria",       "Niamey",        13.52),
        ("Senegal",       "Dakar",         14.69),
        ("Sierra Leone",  "Freetown",       8.48),
        ("South Africa",  "Capetown",     -33.93),
        ("Togo",          "Lome",           6.17),
        ("Tunisia",       "Tunis",         36.80),
        ("Tanzania",      "Dar Es Salaam", -6.79),
        ("Egypt",         "Cairo",         30.04),
        ("Kenya",         "Nairobi",       -1.29),
        ("Ghana",         "Accra",          5.55),
        ("Ethiopia",      "Addis Ababa",    8.98),
        ("South Africa",  "Johannesburg", -26.20),
    ],
    "North America": [
        ("United States", "New York",      40.71),
        ("United States", "Los Angeles",   34.05),
        ("United States", "Chicago",       41.88),
        ("Canada",        "Toronto",       43.65),
        ("Mexico",        "Mexico City",   19.43),
    ],
    "South America": [
        ("Brazil",        "São Paulo",    -23.55),
        ("Argentina",     "Buenos Aires", -34.60),
        ("Chile",         "Santiago",     -33.45),
        ("Peru",          "Lima",         -12.04),
        ("Colombia",      "Bogotá",         4.71),
    ],
    "Oceania": [
        ("Australia",     "Sydney",       -33.87),
        ("Australia",     "Melbourne",    -37.81),
        ("New Zealand",   "Auckland",     -36.85),
    ],
}
size_days = {
    "min": 365*10,
    "small": 365 * 20,
    "full": 365 * 30,
}


def synthetic_temperature(date: datetime, latitude: float) -> float:
    day_of_year = date.timetuple().tm_yday
    peak_day = 196 if latitude >= 0 else 15
    baseline = 25 - 0.30 * abs(latitude)
    amplitude = 7 + 0.10 * abs(latitude)
    seasonal = amplitude * cos(2 * pi * (day_of_year - peak_day) / 365)
    noise = random.gauss(0, 2)
    return round(baseline + seasonal + noise, 1)


def generate_dataset(start_date: datetime, num_days: int, out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8") as fh:
        for continent, locs in continent_locations.items():
            for day_offset in range(num_days):
                date = start_date + timedelta(days=day_offset)
                for country, city, lat in locs:
                    temp = synthetic_temperature(date, lat)
                    fh.write(
                        f"{date.month}\t{date.day}\t{date.year}\t"
                        f"{temp}\t{country}\t{city}\t{continent}\n"
                    )


if __name__ == "__main__":
    random.seed(42)
    parser = argparse.ArgumentParser(
        description="Generate consecutive, continent-grouped temperature data."
    )
    parser.add_argument("output_file", help="Destination .txt file")
    parser.add_argument(
        "--size",
        choices=size_days,
        default="min",
        help="Dataset scale: min | small | mid | full",
    )
    parser.add_argument(
        "--start_date", default="1995-01-01", help="Start date (YYYY-MM-DD)"
    )

    args = parser.parse_args()
    start_dt = datetime.strptime(args.start_date, "%Y-%m-%d")
    num_days = size_days[args.size]

    generate_dataset(start_dt, num_days, args.output_file)

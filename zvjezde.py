from datetime import datetime,timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import AltAz,EarthLocation,Galactic,SkyCoord
from astropy.time import Time
from matplotlib.patches import Circle,Rectangle
from skyfield.api import Star,load,wgs84
from skyfield.data import hipparcos

#=========================
#USER SETTINGS
#=========================
MODE="gift"                  #"gift" or "exact"
ENVIRONMENT="dark"           #"city","suburban","dark"

TITLE="EXAMPLE TITLE"
CITY="Your_city"
LAT=Your_latitude
LON=Your_longitude

#Enter your desired date and time
DT_LOCAL=datetime(2026,2,28,22,49,0,tzinfo=ZoneInfo("Europe/Sarajevo"))

SHOW_CONSTELLATIONS=True
CONSTELLATION_FILE="constellationship.fab"

#=========================
#PRINT SETTINGS
#=========================
PRINT_FORMAT="A3"            #"A4" or "A3"
PRINT_DPI=1200                #300,600 or 1200

EXPORT_PNG=True
EXPORT_PDF=False
EXPORT_TIFF=False

BASE_FILENAME=f"star_poster_{MODE}_{PRINT_FORMAT}_{PRINT_DPI}dpi"

OUTPUT_DIR=Path(__file__).resolve().parent/"exports"
OUTPUT_DIR.mkdir(parents=True,exist_ok=True)

#=========================
#STYLE
#=========================
def get_style(mode,environment):
    """
    Return style settings for the poster.
    """
    env_factor={
        "city":0.28,
        "suburban":0.60,
        "dark":1.00,
    }.get(environment,0.28)

    if mode=="gift":
        return{
            "star_limit":6.2,
            "constellation_alpha":0.50,
            "constellation_lw":0.35,
            "mw_alpha_scale":1.00,
            "mw_alpha_cap":0.16,
            "mw_size_mean":11.0,
            "mw_size_std":3.4,
            "mw_size_min":2.0,
            "mw_size_max":26.0,
            "bright_glow":True,
            "bright_glow_alpha":0.10,
            "bright_glow_scale":5.5,
            "rim_1_alpha":0.055,
            "rim_2_alpha":0.035,
            "star_size_max":34.0,
            "star_color_mode":"gift",
        }

    return{
        "star_limit":4.4 if environment=="city" else 5.4 if environment=="suburban" else 6.0,
        "constellation_alpha":0.22,
        "constellation_lw":0.30,
        "mw_alpha_scale":0.45*env_factor,
        "mw_alpha_cap":0.035 if environment=="city" else 0.07 if environment=="suburban" else 0.12,
        "mw_size_mean":7.0,
        "mw_size_std":2.0,
        "mw_size_min":1.2,
        "mw_size_max":10.0,
        "bright_glow":False,
        "bright_glow_alpha":0.04,
        "bright_glow_scale":2.5,
        "rim_1_alpha":0.025,
        "rim_2_alpha":0.015,
        "star_size_max":20.0,
        "star_color_mode":"exact",
    }


STYLE=get_style(MODE,ENVIRONMENT)
STAR_LIMIT=STYLE["star_limit"]

#=========================
#HELPERS
#=========================
def wrapdiff_deg(value,center):
    """
    Smallest wrapped angular difference in degrees.
    """
    return ((value-center+180)%360)-180


def project_altaz(alt_deg,az_deg):
    """
    Convert altitude/azimuth into 2D coordinates inside the sky circle.
    """
    radius=(90-alt_deg)/90.0
    theta=np.deg2rad(az_deg)
    x_coord=radius*np.sin(theta)
    y_coord=radius*np.cos(theta)
    return x_coord,y_coord


def project_galactic_visible(l_vals,b_vals,location,obstime):
    """
    Convert Galactic coordinates to local AltAz and keep only visible points.
    Returns None if nothing is above the horizon.
    """
    coords=SkyCoord(l=l_vals*u.deg,b=b_vals*u.deg,frame=Galactic())
    altaz_coords=coords.transform_to(AltAz(obstime=obstime,location=location))
    visible_mask=altaz_coords.alt.deg>0

    if not np.any(visible_mask):
        return None

    alt_visible=altaz_coords.alt.deg[visible_mask]
    az_visible=altaz_coords.az.deg[visible_mask]
    x_coord,y_coord=project_altaz(alt_visible,az_visible)

    return{
        "x":x_coord,
        "y":y_coord,
        "mask":visible_mask,
        "alt":alt_visible,
        "az":az_visible,
    }


def deg_to_dms_text(lat,lon):
    """
    Convert decimal coordinates into degrees/minutes/seconds text.
    """
    def format_one(value,is_lat=True):
        hemisphere=("N" if value>=0 else "S") if is_lat else ("E" if value>=0 else "W")
        abs_value=abs(value)
        degrees=int(abs_value)
        minutes_float=(abs_value-degrees)*60
        minutes=int(minutes_float)
        seconds=int(round((minutes_float-minutes)*60))

        if seconds==60:
            seconds=0
            minutes+=1
        if minutes==60:
            minutes=0
            degrees+=1

        return f'{degrees}° {minutes}\' {seconds}" {hemisphere}'

    return f"{format_one(lat,True)} {format_one(lon,False)}"


def get_print_size(fmt):
    """
    Return physical page size in inches for Matplotlib export.
    """
    fmt=fmt.upper()
    if fmt=="A4":
        return 8.27,11.69
    if fmt=="A3":
        return 11.69,16.54
    raise ValueError("PRINT_FORMAT mora biti 'A4' ili 'A3'")


def save_outputs(fig,output_dir,base_filename,print_dpi):
    """
    Save poster in selected output formats.
    """
    png_file=output_dir/f"{base_filename}.png"
    pdf_file=output_dir/f"{base_filename}.pdf"
    tiff_file=output_dir/f"{base_filename}.tiff"

    if EXPORT_PNG:
        fig.savefig(
            str(png_file),
            dpi=print_dpi,
            facecolor="black",
            bbox_inches=None,
            pad_inches=0,
        )

    if EXPORT_PDF:
        fig.savefig(
            str(pdf_file),
            dpi=print_dpi,
            facecolor="black",
            bbox_inches=None,
            pad_inches=0,
        )

    if EXPORT_TIFF:
        fig.savefig(
            str(tiff_file),
            dpi=print_dpi,
            facecolor="black",
            bbox_inches=None,
            pad_inches=0,
        )

    return png_file,pdf_file,tiff_file


def star_colors_by_mode(visible_df,hip_ids,mag,mode):
    """
    Generate star colors and transparency.
    """
    if mode=="gift":
        palette=np.array([
            [0.96,0.98,1.00],
            [1.00,0.97,0.90],
            [0.82,0.89,1.00],
            [1.00,0.90,0.62],
            [0.92,0.94,1.00],
        ])
        color_idx=np.mod(hip_ids.astype(int),len(palette))
        colors=palette[color_idx].copy()

        faint_factor=np.clip((mag-1.5)/4.5,0,1)
        colors=colors*(1-faint_factor[:,None]*0.60)+1.0*(faint_factor[:,None]*0.60)

        alpha=np.clip(1.0-(mag/8.0),0.38,0.96)
        return np.column_stack([colors,alpha])

    if "bv_magnitude" in visible_df.columns:
        bv=visible_df["bv_magnitude"].to_numpy()
        bv=np.nan_to_num(bv,nan=0.65)
    else:
        bv=np.full(len(mag),0.65)

    mix=np.clip((bv+0.3)/2.2,0,1)
    cool=np.array([0.92,0.96,1.00])
    warm=np.array([1.00,0.95,0.88])

    rgb=(1-mix[:,None])*cool+mix[:,None]*warm
    rgb=0.72*1.0+0.28*rgb

    alpha=np.clip(0.92-(mag/8.5),0.22,0.88)
    return np.column_stack([rgb,alpha])


def draw_milky_way(ax,sky_circle,lat,lon,dt_local,style,mode):
    """
    Render the Milky Way using a layered particle approach.
    """
    rng=np.random.default_rng(20260403)

    location=EarthLocation(lat=lat*u.deg,lon=lon*u.deg,height=0*u.m)
    obstime=Time(dt_local.astimezone(timezone.utc))

    if mode=="gift":
        n_wide=90000
        n_mid=70000
        n_core=32000
        wide_sigma=10.5
        mid_sigma=4.8
        core_l_sigma=20.0
        core_b_sigma=7.0
    else:
        n_wide=42000
        n_mid=28000
        n_core=10000
        wide_sigma=8.0
        mid_sigma=3.4
        core_l_sigma=14.0
        core_b_sigma=5.0

    l_wide=rng.uniform(0,360,n_wide)
    b_wide=rng.normal(0,wide_sigma,n_wide)

    l_mid=rng.uniform(0,360,n_mid)
    b_mid=rng.normal(0,mid_sigma,n_mid)

    l_core=(rng.normal(0,core_l_sigma,n_core)+360)%360
    b_core=rng.normal(0,core_b_sigma,n_core)

    gal_l=np.concatenate([l_wide,l_mid,l_core])
    gal_b=np.concatenate([b_wide,b_mid,b_core])

    projected=project_galactic_visible(gal_l,gal_b,location,obstime)
    if projected is None:
        return

    x_coord=projected["x"]
    y_coord=projected["y"]
    visible_mask=projected["mask"]

    gal_l_visible=gal_l[visible_mask]
    gal_b_visible=gal_b[visible_mask]

    density_wide=np.exp(-(np.abs(gal_b_visible)/9.0)**2)
    density_mid=np.exp(-(np.abs(gal_b_visible)/4.2)**2)

    galactic_center_dist=np.abs(wrapdiff_deg(gal_l_visible,0))
    bulge=np.exp(-(galactic_center_dist/18.0)**2)*np.exp(-(np.abs(gal_b_visible)/7.5)**2)

    phase_1=(np.sin(np.deg2rad(gal_l_visible*1.7))+1)/2
    phase_2=(np.sin(np.deg2rad(gal_l_visible*4.9+40))+1)/2
    phase_3=(np.cos(np.deg2rad(gal_l_visible*9.3-15))+1)/2

    texture=0.50+0.28*phase_1+0.14*phase_2+0.08*phase_3
    texture=np.clip(texture,0.35,1.15)

    intensity=(0.30*density_wide+0.85*density_mid+1.15*bulge)*texture
    alpha=np.clip(0.010+style["mw_alpha_scale"]*0.125*intensity,0.003,style["mw_alpha_cap"])

    size=np.clip(
        rng.lognormal(mean=np.log(style["mw_size_mean"]),sigma=0.34,size=len(x_coord)),
        style["mw_size_min"],
        style["mw_size_max"],
    )

    cool=np.array([0.52,0.66,0.98])
    neutral=np.array([0.90,0.92,1.00])
    warm=np.array([1.00,0.80,0.58])

    warm_core=np.clip(0.88*bulge+0.08*phase_1,0,1)
    blue_halo=np.clip(0.75*density_wide*(1-bulge)+0.12*phase_2,0,1)

    rgb=(
        0.42*neutral
        +(0.36*blue_halo)[:,None]*cool
        +(0.44*warm_core)[:,None]*warm
    )
    rgb=np.clip(rgb,0,1)

    glow=ax.scatter(
        x_coord,
        y_coord,
        s=size,
        c=np.column_stack([rgb,alpha]),
        linewidths=0,
        zorder=1,
        rasterized=True,
    )
    glow.set_clip_path(sky_circle)

    n_spark=24000 if mode=="gift" else 7000
    spark_l=rng.uniform(0,360,n_spark)
    spark_b=rng.normal(0,5.0 if mode=="gift" else 3.0,n_spark)

    projected_sparks=project_galactic_visible(spark_l,spark_b,location,obstime)
    if projected_sparks is None:
        return

    spark_x=projected_sparks["x"]
    spark_y=projected_sparks["y"]
    spark_mask=projected_sparks["mask"]

    spark_l_visible=spark_l[spark_mask]
    spark_b_visible=spark_b[spark_mask]

    gc=np.exp(-(np.abs(wrapdiff_deg(spark_l_visible,0))/22.0)**2)*np.exp(-(np.abs(spark_b_visible)/6.5)**2)
    density=np.exp(-(np.abs(spark_b_visible)/3.8)**2)

    spark_alpha=np.clip(0.010+0.030*density+0.050*gc,0.008,0.08 if mode=="gift" else 0.035)
    spark_size=np.clip(rng.lognormal(mean=1.10,sigma=0.38,size=len(spark_x)),1.2,9.0)

    spark_rgba=np.column_stack([
        np.full(len(spark_x),0.96),
        np.full(len(spark_x),0.94),
        np.full(len(spark_x),1.00),
        spark_alpha,
    ])

    spark_plot=ax.scatter(
        spark_x,
        spark_y,
        s=spark_size,
        c=spark_rgba,
        linewidths=0,
        zorder=1.15,
        rasterized=True,
    )
    spark_plot.set_clip_path(sky_circle)


def draw_constellations(ax,sky_circle,hip_ids,x_coord,y_coord,style):
    """
    Draw constellation connection lines from Stellarium data.
    """
    if not SHOW_CONSTELLATIONS:
        return

    constellation_path=Path(__file__).resolve().parent/CONSTELLATION_FILE
    if not constellation_path.exists():
        print(f"[INFO]Fajl '{constellation_path}' nije pronadjen.Preskacem sazvijezdja.")
        return

    from skyfield.data import stellarium

    with open(constellation_path,"rb") as constellation_file:
        constellations=stellarium.parse_constellations(constellation_file)

    hip_to_xy={int(hip):(x_val,y_val) for hip,x_val,y_val in zip(hip_ids,x_coord,y_coord)}

    for _,edges in constellations:
        for hip_a,hip_b in edges:
            point_a=hip_to_xy.get(hip_a)
            point_b=hip_to_xy.get(hip_b)
            if point_a and point_b:
                line=ax.plot(
                    [point_a[0],point_b[0]],
                    [point_a[1],point_b[1]],
                    color=(1,1,1,style["constellation_alpha"]),
                    linewidth=style["constellation_lw"],
                    zorder=3,
                )[0]
                line.set_clip_path(sky_circle)


def draw_nebulae_and_dust(ax,sky_circle,lat,lon,dt_local,mode):
    """
    Add artistic nebula clouds and dark dust lanes along the Milky Way.
    """
    rng=np.random.default_rng(20260404)

    location=EarthLocation(lat=lat*u.deg,lon=lon*u.deg,height=0*u.m)
    obstime=Time(dt_local.astimezone(timezone.utc))

    if mode=="gift":
        centers_l=np.array([350,8,20,34,48,74,96,118,210,248,308,334])
        centers_b=np.array([1,-1,2,-2,1,2,-2,1,1,-2,1,0])
        points_per_cloud=1400
        cloud_sigma_l=4.8
        cloud_sigma_b=1.8
        alpha_base=0.020
    else:
        centers_l=np.array([350,18,42,96,248,334])
        centers_b=np.array([0,-1,1,-1,1,0])
        points_per_cloud=450
        cloud_sigma_l=2.8
        cloud_sigma_b=1.0
        alpha_base=0.010

    nebula_palette=np.array([
        [1.00,0.58,0.66],
        [0.62,0.78,1.00],
        [0.84,0.62,1.00],
        [1.00,0.74,0.46],
        [0.70,0.88,1.00],
    ])

    for index,(center_l,center_b) in enumerate(zip(centers_l,centers_b)):
        cloud_l=(rng.normal(center_l,cloud_sigma_l,points_per_cloud)+360)%360
        cloud_b=rng.normal(center_b,cloud_sigma_b,points_per_cloud)

        projected=project_galactic_visible(cloud_l,cloud_b,location,obstime)
        if projected is None:
            continue

        x_coord=projected["x"]
        y_coord=projected["y"]
        visible_mask=projected["mask"]

        l_visible=cloud_l[visible_mask]
        b_visible=cloud_b[visible_mask]

        dist=np.abs(wrapdiff_deg(l_visible,center_l))
        strength=np.exp(-(dist/(cloud_sigma_l*1.15))**2)*np.exp(-(np.abs(b_visible-center_b)/(cloud_sigma_b*1.8))**2)

        color=nebula_palette[index%len(nebula_palette)]
        size=np.clip(
            rng.lognormal(mean=3.00 if mode=="gift" else 2.55,sigma=0.42,size=len(x_coord)),
            8,
            120 if mode=="gift" else 50,
        )
        alpha=np.clip(alpha_base*(0.35+1.20*strength),0.002,0.045 if mode=="gift" else 0.018)

        rgba=np.column_stack([
            np.full(len(x_coord),color[0]),
            np.full(len(x_coord),color[1]),
            np.full(len(x_coord),color[2]),
            alpha,
        ])

        cloud=ax.scatter(
            x_coord,
            y_coord,
            s=size,
            c=rgba,
            linewidths=0,
            zorder=2.0,
            rasterized=True,
        )
        cloud.set_clip_path(sky_circle)

    dust_centers_l=np.array([352,6,16,30,44,84,104,248,330])
    dust_centers_b=np.array([0,1,-1,1,-1,0,-1,1,0])

    for center_l,center_b in zip(dust_centers_l,dust_centers_b):
        n_points=2600 if mode=="gift" else 900
        dust_l=(rng.normal(center_l,5.4 if mode=="gift" else 3.0,n_points)+360)%360
        dust_b=rng.normal(center_b,0.95 if mode=="gift" else 0.55,n_points)

        projected=project_galactic_visible(dust_l,dust_b,location,obstime)
        if projected is None:
            continue

        x_coord=projected["x"]
        y_coord=projected["y"]

        size=np.clip(
            rng.lognormal(mean=3.18 if mode=="gift" else 2.68,sigma=0.36,size=len(x_coord)),
            10,
            155 if mode=="gift" else 65,
        )
        alpha=np.clip(
            rng.normal(0.040 if mode=="gift" else 0.018,0.006,size=len(x_coord)),
            0.010,
            0.060 if mode=="gift" else 0.028,
        )

        rgba=np.column_stack([
            np.full(len(x_coord),0.018),
            np.full(len(x_coord),0.020),
            np.full(len(x_coord),0.030),
            alpha,
        ])

        dust=ax.scatter(
            x_coord,
            y_coord,
            s=size,
            c=rgba,
            linewidths=0,
            zorder=2.15,
            rasterized=True,
        )
        dust.set_clip_path(sky_circle)

    if mode=="gift":
        lane_l=(rng.normal(0,20,6000)+360)%360
        lane_b=rng.normal(0,1.0,6000)

        projected=project_galactic_visible(lane_l,lane_b,location,obstime)
        if projected is not None:
            x_coord=projected["x"]
            y_coord=projected["y"]

            size=np.clip(rng.lognormal(mean=2.8,sigma=0.38,size=len(x_coord)),8,80)
            alpha=np.clip(rng.normal(0.026,0.005,size=len(x_coord)),0.010,0.040)

            rgba=np.column_stack([
                np.full(len(x_coord),0.025),
                np.full(len(x_coord),0.022),
                np.full(len(x_coord),0.030),
                alpha,
            ])

            lane=ax.scatter(
                x_coord,
                y_coord,
                s=size,
                c=rgba,
                linewidths=0,
                zorder=2.2,
                rasterized=True,
            )
            lane.set_clip_path(sky_circle)

#=========================
#MAIN
#=========================
def make_poster():
    """
    Main function:
    1. Load star data
    2. Calculate visible sky for chosen place and time
    3. Draw poster layers
    4. Export final artwork
    """
    ts=load.timescale()
    t=ts.from_datetime(DT_LOCAL.astimezone(timezone.utc))

    eph=load("de421.bsp")
    earth=eph["earth"]
    observer=wgs84.latlon(LAT,LON)

    with load.open(hipparcos.URL) as catalog_file:
        df=hipparcos.load_dataframe(catalog_file)

    df=df[df["ra_degrees"].notnull()]
    df=df[df["dec_degrees"].notnull()]
    df=df[df["magnitude"].notnull()]
    df=df[df["magnitude"]<=STAR_LIMIT]

    stars=Star.from_dataframe(df)
    apparent=(earth+observer).at(t).observe(stars).apparent()
    alt,az,_=apparent.altaz()

    alt_deg=alt.degrees
    az_deg=az.degrees

    visible_mask=alt_deg>0
    visible_df=df[visible_mask].copy()

    alt_visible=alt_deg[visible_mask]
    az_visible=az_deg[visible_mask]
    x_coord,y_coord=project_altaz(alt_visible,az_visible)

    magnitudes=visible_df["magnitude"].to_numpy()
    hip_ids=visible_df.index.to_numpy()

    star_sizes=np.clip((7.2-magnitudes)**2,1.2,STYLE["star_size_max"])
    star_colors=star_colors_by_mode(visible_df,hip_ids,magnitudes,STYLE["star_color_mode"])

    fig_width,fig_height=get_print_size(PRINT_FORMAT)
    fig=plt.figure(figsize=(fig_width,fig_height),facecolor="black")
    ax=fig.add_axes((0.08,0.34,0.84,0.58))
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    ax.set_xlim(-1.08,1.08)
    ax.set_ylim(-1.08,1.08)
    ax.axis("off")

    sky_circle=Circle((0,0),1,fill=False,edgecolor="white",linewidth=1.4,zorder=10)
    ax.add_patch(sky_circle)

    for radius,line_width,alpha in ((0.985,18,STYLE["rim_1_alpha"]),(0.955,10,STYLE["rim_2_alpha"])):
        rim=Circle((0,0),radius,fill=False,edgecolor=(1,1,1,alpha),linewidth=line_width,zorder=2)
        ax.add_patch(rim)

    draw_milky_way(ax,sky_circle,LAT,LON,DT_LOCAL,STYLE,MODE)
    draw_nebulae_and_dust(ax,sky_circle,LAT,LON,DT_LOCAL,MODE)

    stars_plot=ax.scatter(
        x_coord,
        y_coord,
        s=star_sizes,
        c=star_colors,
        linewidths=0,
        zorder=4,
        rasterized=True,
    )
    stars_plot.set_clip_path(sky_circle)

    if STYLE["bright_glow"]:
        bright_mask=magnitudes<1.8
        if bright_mask.any():
            glow_rgba=np.column_stack([
                np.ones(bright_mask.sum()),
                np.ones(bright_mask.sum()),
                np.ones(bright_mask.sum()),
                np.full(bright_mask.sum(),STYLE["bright_glow_alpha"]),
            ])
            bright_plot=ax.scatter(
                x_coord[bright_mask],
                y_coord[bright_mask],
                s=star_sizes[bright_mask]*STYLE["bright_glow_scale"],
                c=glow_rgba,
                linewidths=0,
                zorder=3,
                rasterized=True,
            )
            bright_plot.set_clip_path(sky_circle)

    draw_constellations(ax,sky_circle,hip_ids,x_coord,y_coord,STYLE)

    ax.add_patch(Circle((0,0),1,fill=False,edgecolor="white",linewidth=1.5,zorder=12))

    fig.add_artist(Rectangle(
        (0.05,0.04),0.90,0.92,
        transform=fig.transFigure,
        fill=False,
        edgecolor="white",
        linewidth=1.6,
    ))

    coord_text=deg_to_dms_text(LAT,LON)
    date_text=DT_LOCAL.strftime("%B %d, %Y").upper()

    fig.text(
        0.5,0.27,TITLE,
        color="white",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        family="serif",
    )

    fig.text(
        0.5,0.17,CITY.upper(),
        color="white",
        ha="center",
        fontsize=13,
        fontweight="bold",
        family="serif",
    )

    fig.text(
        0.5,0.145,coord_text,
        color="white",
        ha="center",
        fontsize=11,
        fontweight="bold",
        family="serif",
    )

    fig.text(
        0.5,0.12,date_text,
        color="white",
        ha="center",
        fontsize=11,
        fontweight="bold",
        family="serif",
    )

    png_file,pdf_file,tiff_file=save_outputs(fig,OUTPUT_DIR,BASE_FILENAME,PRINT_DPI)

    plt.show()

    print("Poster saved as:")
    if EXPORT_PNG:
        print(f"- {png_file}")
    if EXPORT_PDF:
        print(f"- {pdf_file}")
    if EXPORT_TIFF:
        print(f"- {tiff_file}")


if __name__=="__main__":
    make_poster()
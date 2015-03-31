#include "media-io\format-conversion.h"

#include <obs-module.h>
#include <util/circlebuf.h>

#include <algorithm>
#include <cstdint>
#include <immintrin.h>

#ifndef SEC_TO_NSEC
#define SEC_TO_NSEC 1000000000ULL
#endif

#ifndef MSEC_TO_NSEC
#define MSEC_TO_NSEC 1000000ULL
#endif

#define SIZE_OSF sizeof(struct obs_source_frame*)

struct deinterlace_data {
	obs_source_t                   *context;

	/* contains struct obs_source_frame* */
	struct circlebuf               video_frames;

  struct obs_source_frame output_frame;

  bool reset_video;
};

static const char *deinterlace_filter_name(void)
{
	return obs_module_text("DeinterlaceFilter");
}

static void free_video_data(struct deinterlace_data *filter,
		obs_source_t *parent)
{
	while (filter->video_frames.size) {
		struct obs_source_frame *frame;

		circlebuf_pop_front(&filter->video_frames, &frame,
				sizeof(struct obs_source_frame*));
		obs_source_release_frame(parent, frame);
	}
}

static inline void free_audio_packet(struct obs_audio_data *audio)
{
	for (size_t i = 0; i < MAX_AV_PLANES; i++)
		bfree(audio->data[i]);
	memset(audio, 0, sizeof(*audio));
}

static void deinterlace_filter_update(void *data, obs_data_t *settings)
{
	deinterlace_data *filter = (deinterlace_data *)data;
  filter->reset_video = true;
}

static void *deinterlace_filter_create(obs_data_t *settings,
		obs_source_t *context)
{
	deinterlace_data *filter = (deinterlace_data *)bzalloc(sizeof(*filter));

	filter->context = context;
  circlebuf_reserve(&filter->video_frames, 3 * sizeof(struct obs_source_frame*));
	deinterlace_filter_update(filter, settings);

	return filter;
}

static void deinterlace_filter_destroy(void *data)
{
	deinterlace_data *filter = (deinterlace_data *)data;

	circlebuf_free(&filter->video_frames);
	bfree(data);
}

static obs_properties_t *deinterlace_filter_properties(void *data)
{
	obs_properties_t *props = obs_properties_create();
	UNUSED_PARAMETER(data);
	return props;
}

static void deinterlace_filter_remove(void *data, obs_source_t *parent)
{
	deinterlace_data *filter = (deinterlace_data *)data;

	free_video_data(filter, parent);
}

static void i420_to_rgba(struct obs_source_frame *dst,
		const struct obs_source_frame *src) {

}

static void deinterlace_plane(uint8_t *dst, const uint8_t *src, uint32_t linesize, uint32_t height) {
  for (int i = 0; i < height; i += 2) {
    memcpy(dst + linesize * i, src + linesize * i, linesize);
    memcpy(dst + linesize * (i + 1), src + linesize * i, linesize);
  }
}

// =========== GNU General Public License code start =================

//#define MIN(a,b) ((a) > (b) ? (b) : (a))
//#define MAX(a,b) ((a) < (b) ? (b) : (a))
//#define ABS(a) ((a) > 0 ? (a) : (-(a)))
#define MIN(a,b) std::min(a,b)
#define MAX(a,b) std::max(a,b)
#define ABS(a) std::abs(a)

#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

inline float halven(float f) { return f*0.5f; }
inline int halven(int i) { return i>>1; }

inline int one1(unsigned char *) { return 1; }
inline float one1(float *) { return 1.0f/256.0f; }

template<int ch,typename Comp,typename Diff>
inline void filter_line_c(int mode, Comp *dst, const Comp *prev, const Comp *cur, const Comp *next, int width, long pitch, int parity){
    const Comp *prev2 = parity ? prev : cur ;
    const Comp *next2 = parity ? cur  : next;

    for(int x = 0; x < width; ++x){
        Diff c_above = cur[-pitch];
        Diff d = halven(prev2[0] + next2[0]);
        Diff c_below = cur[+pitch];
        Diff temporal_diff0=halven( ABS(prev2[0] - next2[0]) );
        Diff temporal_diff1=halven( ABS(prev[-pitch] - c_above) + ABS(prev[+pitch] - c_below) );
        Diff temporal_diff2=halven( ABS(next[-pitch] - c_above) + ABS(next[+pitch] - c_below) );
        Diff diff= MAX3(temporal_diff0, temporal_diff1, temporal_diff2);
        Diff spatial_pred= halven(c_above + c_below);
        Diff spatial_score = ABS(cur[-pitch - ch] - cur[+pitch - ch])
                           + ABS(c_above - c_below)
                           + ABS(cur[-pitch + ch] - cur[+pitch + ch])
                           - one1((Comp*)0);

#define CHECK(j)\
    {   Diff score= ABS(cur[-pitch-ch+ j] - cur[+pitch-ch- j])\
                 + ABS(cur[-pitch   + j] - cur[+pitch   - j])\
                 + ABS(cur[-pitch+ch+ j] - cur[+pitch+ch- j]);\
        if(score < spatial_score){\
            spatial_score= score;\
            spatial_pred= halven(cur[-pitch  + j] + cur[+pitch  - j]);\

        CHECK(-ch) CHECK(-(ch*2)) }} }}
        CHECK( ch) CHECK( (ch*2)) }} }}

        if(mode<2){
            Diff b= halven(prev2[-2*pitch] + next2[-2*pitch]);
            Diff f= halven(prev2[+2*pitch] + next2[+2*pitch]);

            Diff max= MAX3(d-c_below, d-c_above, MIN(b-c_above, f-c_below));
            Diff min= MIN3(d-c_below, d-c_above, MAX(b-c_above, f-c_below));

            diff= MAX3(diff, min, -max);
        }

        if(spatial_pred > d + diff)
           spatial_pred = d + diff;
        else if(spatial_pred < d - diff)
           spatial_pred = d - diff;
        
        dst[0] = (Comp)spatial_pred;

        dst+=ch;
        cur+=ch;
        prev+=ch;
        next+=ch;
        prev2+=ch;
        next2+=ch;
    }
}

typedef unsigned char c16 __attribute__((ext_vector_type(16)));
typedef short s16 __attribute__((ext_vector_type(16)));

static s16 loads(const std::uint8_t *__restrict addr) {
  return __builtin_convertvector(*(const c16 *)addr, s16);
}

static s16 abs_diff(s16 a, s16 b) {
  return _mm256_abs_epi16(a - b);
}

static s16 avg_of_abs_diff(s16 a, s16 b) {
  return abs_diff(a, b) >> 1; 
}

static s16 calc_spatial_score(const std::uint8_t *cur, int pitch, int j) {
  s16 cur_above_behind = loads(cur - pitch - 1 + j);
  s16 cur_above = loads(cur - pitch + j);
  s16 cur_above_forward = loads(cur - pitch + 1 + j);
  s16 cur_below_behind = loads(cur + pitch - 1 - j);
  s16 cur_below = loads(cur + pitch - j);
  s16 cur_below_forward = loads(cur + pitch + 1 - j);
  return abs_diff(cur_above_behind, cur_below_behind) +
         abs_diff(cur_above, cur_below) +
         abs_diff(cur_above_forward, cur_below_forward);
}

static s16 clamp(s16 val, s16 max, s16 min) {
  return _mm256_min_epi16(_mm256_max_epi16(val, min), max);
}

static s16 calc_spatial_pred(const std::uint8_t *cur, int pitch, int j) {
  return _mm256_avg_epu16(loads(cur - pitch + j), loads(cur + pitch - j));
}

static void filter_line_c(int mode, std::uint8_t *__restrict dst, const std::uint8_t *__restrict prev, const std::uint8_t *__restrict cur, const std::uint8_t *__restrict next, int width, long pitch, int parity) {
    const std::uint8_t *prev2 = parity ? prev : cur ;
    const std::uint8_t *next2 = parity ? cur  : next;

    for (int x = 0; x < width; x += sizeof(__m128i)) {
      s16 c_above = loads(cur - pitch);
      s16 c_below = loads(cur + pitch);
      s16 c_before =loads(prev2);
      s16 c_after = loads(next2);
      s16 d = _mm256_avg_epu16(c_before, c_after);
      s16 temporal_diff0 = avg_of_abs_diff(c_before, c_after);
      s16 c_prev_above = loads(prev + pitch);
      s16 c_prev_below = loads(prev - pitch);
      s16 c_next_above = loads(next + pitch);
      s16 c_next_below = loads(next - pitch);
      s16 temporal_diff1 = _mm256_avg_epu16(abs_diff(c_prev_above, c_above), abs_diff(c_prev_below, c_below));
      s16 temporal_diff2 = _mm256_avg_epu16(abs_diff(c_next_above, c_above), abs_diff(c_next_below, c_below));
      s16 diff = _mm256_max_epi16(_mm256_max_epi16(temporal_diff0, temporal_diff1), temporal_diff2);
      s16 spatial_pred = _mm256_avg_epu16(c_above, c_below);
      s16 spatial_score = calc_spatial_score(cur, pitch, 0) - 1;
      
#define CHECK(j)\
      {\
        s16 score = calc_spatial_score(cur, pitch, j);\
        s16 new_spatial_pred = calc_spatial_pred(cur, pitch, j);\
        s16 mask = _mm256_cmpgt_epi16(spatial_score, score);\
        score = _mm256_blendv_epi8(spatial_score, score, mask);\
        spatial_pred = _mm256_blendv_epi8(spatial_pred, new_spatial_pred, mask);\
      }

      CHECK(-1)
      CHECK(-2)
      CHECK(1)
      CHECK(2)

      s16 b = _mm256_avg_epu16(loads(prev2 - 2 * pitch), loads(next2 - 2 * pitch));
      s16 f = _mm256_avg_epu16(loads(prev2 + 2 * pitch), loads(next2 + 2 * pitch));
      
      s16 max = _mm256_max_epu16(_mm256_max_epu16(d - c_below, d - c_above), _mm256_min_epu16(b - c_above, f - c_below));
      s16 min = _mm256_min_epu16(_mm256_min_epu16(d - c_below, d - c_above), _mm256_max_epu16(b - c_above, f - c_below));
      
      diff = _mm256_max_epu16(_mm256_max_epu16(diff, min), -max);
      
      s16 dpdiff = _mm256_add_epi16(d, diff);
      s16 dmdiff = _mm256_sub_epi16(d, diff);
      spatial_pred = clamp(spatial_pred, dpdiff, dmdiff);
      
      _mm_storeu_si128((__m128i *)cur, __builtin_convertvector(spatial_pred, c16));
      
      dst += sizeof(__m128i);
      prev += sizeof(__m128i);
      cur += sizeof(__m128i);
      next += sizeof(__m128i);
      prev2 += sizeof(__m128i);
      next2 += sizeof(__m128i);
    }
}

inline void interpolate(unsigned char *dst, const unsigned char *cur0,  const unsigned char *cur2, int w)
{
    int x;
    for (x=0; x<w; x++) {
        dst[x] = (cur0[x] + cur2[x] + 1)>>1; // simple average
    }
}

inline void interpolate(float *dst, const float *cur0,  const float *cur2, int w)
{
    int x;
    for (x=0; x<w; x++) {
        dst[x] = (cur0[x] + cur2[x] )*0.5f; // simple average
    }
}


template<int ch,typename Comp,typename Diff>
static void filter_plane(int mode, Comp *dst, long dst_stride, const Comp *prev0, const Comp *cur0, const Comp *next0, long refs, int w, int h, int parity, int tff)
{
	int y=0;

    if(((y ^ parity) & 1)) {
        memcpy(dst, cur0 + refs, w*ch*sizeof(Comp));// duplicate 1
    }else{
        memcpy(dst, cur0, w*ch*sizeof(Comp));
    }

    y=1;
    
    if(((y ^ parity) & 1)) {
        interpolate(dst + dst_stride, cur0, cur0 + refs*2, w*ch);   // interpolate 0 and 2
    }else{
        memcpy(dst + dst_stride, cur0 + refs, w*ch*sizeof(Comp)); // copy original
    }

    for(y=2; y<h-2; y++) {
        if(((y ^ parity) & 1)) {
            const Comp *prev= prev0 + y*refs;
            const Comp *cur = cur0 + y*refs;
            const Comp *next= next0 + y*refs;
            Comp *dst2= dst + y*dst_stride;

            for(int k=0;k<ch;k++)
                filter_line_c(mode, dst2+k, prev+k, cur+k, next+k, w, refs, (parity ^ tff));
        }else{
            memcpy(dst + y*dst_stride, cur0 + y*refs, w*ch*sizeof(Comp)); // copy original
        }
    }

    y=h-2;

    if(((y ^ parity) & 1)) {
        interpolate(dst + (h-2)*dst_stride, cur0 + (h-3)*refs, cur0 + (h-1)*refs, w*ch);   // interpolate h-3 and h-1
    }else{
      memcpy(dst + (h-2)*dst_stride, cur0 + (h-2)*refs, w*ch*sizeof(Comp)); // copy original
    }

    y=h-1;
    
    if(((y ^ parity) & 1)){
        memcpy(dst + (h-1)*dst_stride, cur0 + (h-2)*refs, w*ch*sizeof(Comp)); // duplicate h-2
    }else{
        memcpy(dst + (h-1)*dst_stride, cur0 + (h-1)*refs, w*ch*sizeof(Comp)); // copy original
    }
}

// =========== GNU General Public License code end =================

static void deinterlace_frame(obs_source_frame *dst,
    const obs_source_frame *prev, const obs_source_frame *cur,
    const obs_source_frame *next)
{
	dst->flip         = cur->flip;
	dst->full_range   = cur->full_range;
	dst->timestamp    = cur->timestamp;
	memcpy(dst->color_matrix, cur->color_matrix, sizeof(float) * 16);
	if (!dst->full_range) {
		size_t const size = sizeof(float) * 3;
		memcpy(dst->color_range_min, cur->color_range_min, size);
		memcpy(dst->color_range_max, cur->color_range_max, size);
	}

  
  //deinterlace_plane(dst->data[0], cur->data[0], cur->linesize[0], cur->height);
  filter_plane<1, unsigned char, int>(0, dst->data[0], dst->linesize[0], prev->data[0], cur->data[0], next->data[0], cur->linesize[0], cur->width, cur->height, 0, 0);

	switch (cur->format) {
	case VIDEO_FORMAT_I420:

    //deinterlace_plane(dst->data[0], src->data[0], src->linesize[0], src->height);
    deinterlace_plane(dst->data[1], cur->data[1], cur->linesize[1], cur->height / 2);
    deinterlace_plane(dst->data[2], cur->data[2], cur->linesize[2], cur->height / 2);
		break;

	case VIDEO_FORMAT_NV12:
    //deinterlace_plane(dst->data[0], src->data[0], src->linesize[0], src->height);
    deinterlace_plane(dst->data[1], cur->data[1], cur->linesize[1], cur->height / 2);
		break;

	case VIDEO_FORMAT_YVYU:
	case VIDEO_FORMAT_YUY2:
	case VIDEO_FORMAT_UYVY:
	case VIDEO_FORMAT_NONE:
	case VIDEO_FORMAT_RGBA:
	case VIDEO_FORMAT_BGRA:
	case VIDEO_FORMAT_BGRX:
		//deinterlace_plane(dst->data[0], src->data[0], src->linesize[0], src->height);
    break;
	}
}


static struct obs_source_frame *deinterlace_filter_video(void *data,
		struct obs_source_frame *frame)
{
	deinterlace_data *filter = (deinterlace_data *)data;
	obs_source_t *parent = obs_filter_get_parent(filter->context);
  obs_source_frame *latest;
  obs_source_frame *middle;
  obs_source_frame *last;

	if (filter->reset_video) {
    bfree(filter->output_frame.data[0]);
    obs_source_frame_init(&filter->output_frame, frame->format, frame->width, frame->height);
		free_video_data(filter, parent);
		filter->reset_video = false;
	}

  // Don't start deinterlacing until we have 3 frames.
  if (filter->video_frames.size / SIZE_OSF < 3) {
    circlebuf_push_back(&filter->video_frames, &frame, SIZE_OSF);
    return NULL;
  }

  circlebuf_pop_front(&filter->video_frames, NULL, SIZE_OSF);
  circlebuf_push_back(&filter->video_frames, &frame,
			sizeof(struct obs_source_frame*));

  circlebuf_pop_front(&filter->video_frames, &last, SIZE_OSF);
  circlebuf_pop_front(&filter->video_frames, &middle, SIZE_OSF);
  circlebuf_pop_front(&filter->video_frames, &latest, SIZE_OSF);

  deinterlace_frame(&filter->output_frame, last, middle, latest);

  // do magic!
  //for (int i = 0; i < middle->height; i += 2) {
    //memcpy(filter->output_frame.data[0] + frame->linesize[0] * i, middle->data[0] + frame->linesize[0] * i, frame->linesize[0]);
    //memcpy(filter->output_frame.data[0] + frame->linesize[0] * (i + 1), middle->data[0] + frame->linesize[0] * i, frame->linesize[0]);
  //}
  
  // memcpy(filter->output_frame.data[0], middle->data[0], middle->linesize[0] * middle->height);

  circlebuf_push_back(&filter->video_frames, &last, SIZE_OSF);
  circlebuf_push_back(&filter->video_frames, &middle, SIZE_OSF);
  circlebuf_push_back(&filter->video_frames, &latest, SIZE_OSF);

  return &filter->output_frame;
}

static obs_source_info deinterlace_filter;
extern "C" void load_deinterlace_filter() {
	deinterlace_filter.id                            = "deinterlace_filter";
	deinterlace_filter.type                          = OBS_SOURCE_TYPE_FILTER;
	deinterlace_filter.output_flags                  = OBS_SOURCE_VIDEO | OBS_SOURCE_ASYNC;
	deinterlace_filter.get_name                      = deinterlace_filter_name;
	deinterlace_filter.create                        = deinterlace_filter_create;
	deinterlace_filter.destroy                       = deinterlace_filter_destroy;
	deinterlace_filter.update                        = deinterlace_filter_update;
	deinterlace_filter.get_properties                = deinterlace_filter_properties;
	deinterlace_filter.filter_video                  = deinterlace_filter_video;
	deinterlace_filter.filter_remove                 = deinterlace_filter_remove;
  obs_register_source(&deinterlace_filter);
}
